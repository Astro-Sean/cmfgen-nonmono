"""Regression tests for STELLA->CMFGEN preparation consistency.

Covers the hardening changes:
  - strict common SN_HYDRO_DATA reader (grid_mapping.sn_hydro_io)
  - prepare_VADAT.parse_sn_hydro canonical checks + na_cgs propagation
  - proxy_populations_v2 VEL_LAW=11 hydro-reference abundances/ED handling
  - convert_stella_swd_to_cmfgen writer validation, ED source ordering,
    post-write boundary metadata
  - submit script proxy flag plumbing

Run:  pytest scripts/test_cmfgen_preparation_consistency.py -v
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import proxy_populations_v2 as pp
import prepare_VADAT as pv
import convert_stella_swd_to_cmfgen as cvt
from grid_mapping.sn_hydro_io import (
    read_sn_hydro_data, validate_sn_hydro_data)

SCRIPTS = Path(os.path.dirname(os.path.abspath(__file__)))
REPO = SCRIPTS.parent


class TestSobLambdaNewtonProfile:
    def test_default_schedule_and_local_ba_controls(self):
        lines = [
            '2.0D2 [LAM_VAL]\n',
            '4.0D2 [LAM_ENTER]\n',
            '2.0D1 [LAM_EXIT]\n',
            '3 [LAM_N_STABLE]\n',
            '20 [NUM_LAM]\n',
            'T [DO_BA_RES_ACC]\n',
            '1.0D4 [BA_ENTER_RAW]\n',
            '20 [BA_PROBE_LAM]\n',
            '1.3D0 [BA_DAMP_CAP]\n',
        ]
        recommended = pv._build_non_mono_sob_recommended(lines, set())
        generated = pv._apply_forced_vadat_keys(lines, recommended)

        expected = {
            'LAM_VAL': '2.0D2',
            'LAM_ENTER': '2.0D2',
            'LAM_EXIT': '2.0D2',
            'LAM_N_STABLE': '1',
            'NUM_LAM': '2',
            'DO_BA_RES_ACC': 'F',
            'BA_ENTER_RAW': '1.0D30',
            'BA_PROBE_LAM': '0',
            'BA_DAMP_CAP': '0.0D0',
            'FRESH_BA_EXIT': 'T',
            'N_FULL_TREL': '0',
            'COMP_BA': 'T',
            'FIX_BA': '1.0D1',
            'N_FIX_BA': '2',
        }
        for key, value in expected.items():
            assert pv._get_key_token(generated, key) == value

    def test_explicit_ba_controls_are_preserved(self):
        lines = [
            '7.5D1 [LAM_VAL]\n',
            '3.0D2 [LAM_ENTER]\n',
            'T [DO_BA_RES_ACC]\n',
            '5 [NUM_LAM]\n',
        ]
        explicit = {'LAM_VAL', 'LAM_ENTER', 'DO_BA_RES_ACC', 'NUM_LAM'}
        recommended = pv._build_non_mono_sob_recommended(lines, explicit)
        generated = pv._apply_forced_vadat_keys(lines, recommended)

        assert pv._get_key_token(generated, 'LAM_VAL') == '7.5D1'
        assert pv._get_key_token(generated, 'LAM_ENTER') == '3.0D2'
        assert pv._get_key_token(generated, 'LAM_EXIT') == '7.5D1'
        assert pv._get_key_token(generated, 'DO_BA_RES_ACC') == 'T'
        assert pv._get_key_token(generated, 'NUM_LAM') == '5'


def _fmt(v, exp='E'):
    return f'{v:.6E}' if exp == 'E' else f'{v:.6E}'.replace('E', 'D')

def _hydro_lines(nd=4, exp='E', niso=0, shuffle=False, preamble=False,
                 ed=None, kappa=None, fracs=None, r=None):
    r = np.linspace(4.0e4, 1.0e4, nd) if r is None else np.asarray(r, float)
    v = np.linspace(6000.0, 1000.0, nd)
    sigma = np.zeros(nd)
    t = np.linspace(0.6, 2.0, nd)
    rho = np.linspace(1e-12, 1e-9, nd)
    ed = np.linspace(1e8, 1e12, nd) if ed is None else np.asarray(ed, float)
    natom = np.linspace(1e10, 1e14, nd)
    kappa = np.full(nd, 0.34) if kappa is None else np.asarray(kappa, float)
    fracs = {'HYD': np.full(nd, 0.7), 'HE': np.full(nd, 0.3)} if fracs is None else fracs

    def blk(label, arr):
        return [label, '   ' + '   '.join(_fmt(x, exp) for x in arr), '']

    headers = [
        f'Number of data points:{nd:>10d}',
        f'Number of mass fractions:{len(fracs):>10d}',
        f'Number of isotopes:{niso:>10d}',
        f'Time(days) since explosion:     {_fmt(30.0, exp)}',
        f'Energy Lum [lsun]:{_fmt(1.0e6, exp)}',
        f'Inner boundary radius [cm]:{_fmt(1.0e14, exp)}',
        f'Inner boundary temperature [K]:{_fmt(2.0e4, exp)}',
    ]
    if shuffle:
        headers = [headers[i] for i in (3, 0, 5, 1, 6, 2, 4)]
    lines = (['', ''] if preamble else []) + headers + ['']
    lines += blk('Radius grid (10^10cm)', r)
    lines += blk('Velocity (km/s)', v)
    lines += blk('Sigma (dlnV/dlnr-1)', sigma)
    lines += blk('Temperature (10^4 K)', t)
    lines += blk('Density (gm/cm^3)', rho)
    lines += blk('Atom density (/cm^3)', natom)
    lines += blk('Electron density (/cm^3)', ed)
    lines += blk('Kappa (cm^2/gm)', kappa)
    for label, arr in fracs.items():
        lines += blk(f'{label} mass fraction', arr)
    return lines

def _write_hydro(path, **kw):
    Path(path).write_text('\n'.join(_hydro_lines(**kw)) + '\n')
    return str(path)

class TestCanonicalReader:

    def test_shuffled_headers_preamble_d_exponents(self, tmp_path):
        p = _write_hydro(tmp_path / 'SN_HYDRO_DATA', nd=4, exp='D',
                         shuffle=True, preamble=True)
        d = read_sn_hydro_data(p)
        assert d.n_points == 4
        assert d.n_mass_fractions == 2
        assert d.n_isotopes == 0
        assert abs(d.time_days - 30.0) < 1e-6
        np.testing.assert_allclose(d.r_1e10cm, np.linspace(4.0e4, 1.0e4, 4))
        np.testing.assert_allclose(d.electron_density, np.linspace(1e8, 1e12, 4))
        assert abs(float(d.headers['Energy Lum [lsun]'].replace('D', 'E')) - 1e6) < 1.0

    def test_short_block_rejected(self, tmp_path):
        lines = _hydro_lines(nd=4)
        i = lines.index('Radius grid (10^10cm)')
        lines[i + 1] = '   ' + '   '.join(f'{x:.6E}' for x in [4.0e4, 3.0e4, 2.0e4])
        p = tmp_path / 'SN_HYDRO_DATA'
        p.write_text('\n'.join(lines) + '\n')
        with pytest.raises(ValueError):
            read_sn_hydro_data(p)

    def test_surplus_row_rejected(self, tmp_path):
        lines = _hydro_lines(nd=4)
        i = lines.index('Radius grid (10^10cm)')
        lines[i + 1] = ('   ' + '   '.join(f'{x:.6E}' for x in np.linspace(4.0e4, 1.0e4, 4))
                        + '\n   5.000000E+03')
        p = tmp_path / 'SN_HYDRO_DATA'
        p.write_text('\n'.join(lines) + '\n')
        with pytest.raises(ValueError):
            read_sn_hydro_data(p)

    def test_duplicate_radius_block_rejected(self, tmp_path):
        lines = _hydro_lines(nd=4)
        i = lines.index('Velocity (km/s)')
        lines[i:i] = ['Radius grid (10^10cm)',
                      '   ' + '   '.join(f'{x:.6E}' for x in np.linspace(4.0e4, 1.0e4, 4)), '']
        p = tmp_path / 'SN_HYDRO_DATA'
        p.write_text('\n'.join(lines) + '\n')
        with pytest.raises(ValueError, match='Duplicate'):
            read_sn_hydro_data(p)

    def test_duplicate_element_block_rejected(self, tmp_path):
        lines = _hydro_lines(nd=4)
        lines += ['HYD mass fraction',
                  '   ' + '   '.join('7.000000E-01' for _ in range(4)), '']
        p = tmp_path / 'SN_HYDRO_DATA'
        p.write_text('\n'.join(lines) + '\n')
        with pytest.raises(ValueError, match='Duplicate'):
            read_sn_hydro_data(p)

    def test_consumer_rejects_nonmonotonic_radius(self, tmp_path):
        p = _write_hydro(tmp_path / 'SN_HYDRO_DATA', nd=4,
                         r=np.linspace(1.0e4, 4.0e4, 4))
        with pytest.raises(ValueError):
            pp.parse_sn_hydro_data(p)
        with pytest.raises(ValueError):
            pv.parse_sn_hydro(Path(p))

    def test_validator_nonfinite_and_bounds(self, tmp_path):
        p = _write_hydro(tmp_path / 'bad_ed', ed=[1e8, np.nan, 1e10, 1e12])
        rep = validate_sn_hydro_data(p)
        assert not rep.passed and any('Electron density' in e for e in rep.errors)

        p = _write_hydro(tmp_path / 'bad_kappa', kappa=[0.3, -0.1, 0.3, 0.3])
        rep = validate_sn_hydro_data(p)
        assert not rep.passed and any('kappa' in e for e in rep.errors)

        p = _write_hydro(tmp_path / 'bad_frac',
                         fracs={'HYD': np.full(4, 1.2), 'HE': np.full(4, -0.2)})
        rep = validate_sn_hydro_data(p)
        assert not rep.passed and any('[0,1]' in e for e in rep.errors)

class TestRequiredHeaders:

    def test_missing_lum_header_rejected(self, tmp_path):
        lines = [ln for ln in _hydro_lines(nd=4) if not ln.startswith('Energy Lum')]
        p = tmp_path / 'SN_HYDRO_DATA'
        p.write_text('\n'.join(lines) + '\n')
        with pytest.raises(ValueError, match='Energy Lum'):
            pp.parse_sn_hydro_data(str(p))

    def test_explicit_zero_lum_accepted(self, tmp_path):
        lines = _hydro_lines(nd=4)
        i = next(i for i, ln in enumerate(lines) if ln.startswith('Energy Lum'))
        lines[i] = 'Energy Lum [lsun]:   0.000000E+00'
        p = tmp_path / 'SN_HYDRO_DATA'
        p.write_text('\n'.join(lines) + '\n')
        h = pp.parse_sn_hydro_data(str(p))
        assert h.lum_lsun == 0.0

    def test_missing_boundary_headers_rejected(self, tmp_path):
        lines = [ln for ln in _hydro_lines(nd=4) if not ln.startswith('Inner boundary')]
        p = tmp_path / 'SN_HYDRO_DATA'
        p.write_text('\n'.join(lines) + '\n')
        with pytest.raises(ValueError, match='Inner boundary'):
            pp.parse_sn_hydro_data(str(p))


class TestPrepareVadatReader:

    def test_na_cgs_returned_reversed_to_inner_outer(self, tmp_path):
        natom_src = np.geomspace(1e10, 1e14, 4)
        lines = _hydro_lines(nd=4)
        ai = lines.index('Atom density (/cm^3)')
        lines[ai + 1] = '   ' + '   '.join(f'{x:.6E}' for x in natom_src)
        p = tmp_path / 'SN_HYDRO_DATA'
        p.write_text('\n'.join(lines) + '\n')
        hydro = pv.parse_sn_hydro(p)
        assert hydro.na_cgs is not None
        np.testing.assert_allclose(hydro.na_cgs, natom_src[::-1], rtol=1e-5)

    def test_fin_cal_grid_recovers_ed_and_real_na(self, tmp_path):
        nd = 4
        natom_src = np.geomspace(1e10, 1e14, nd)
        ed_src = np.geomspace(1e8, 1e11, nd)
        lines = _hydro_lines(nd=nd, ed=ed_src)
        ai = lines.index('Atom density (/cm^3)')
        lines[ai + 1] = '   ' + '   '.join(f'{x:.6E}' for x in natom_src)
        p = tmp_path / 'SN_HYDRO_DATA'
        p.write_text('\n'.join(lines) + '\n')
        hydro = pv.parse_sn_hydro(p)
        out = tmp_path / 'FIN_CAL_GRID'
        pv.write_fin_cal_grid_from_hydro(hydro, out)
        rows = [ln.split() for ln in out.read_text().splitlines()[1:] if ln.split()]
        assert len(rows) == nd
        cols = np.array([[float(x) for x in row] for row in rows])
        np.testing.assert_allclose(cols[:, 5], natom_src, rtol=2e-4)
        np.testing.assert_allclose(cols[:, 5] * cols[:, 6], ed_src, rtol=2e-4)

    def test_mass_fraction_inline_comments_and_blank_rows(self, tmp_path):
        nd = 4
        he = np.array([0.25, 0.30, 0.30, 0.30])
        lines = _hydro_lines(nd=nd, fracs={'HYD': 1.0 - he, 'HE': he})
        hi = lines.index('HE mass fraction')
        vals = '   ' + '   '.join(f'{x:.6E}' for x in he)
        lines[hi + 1] = '   ' + '   '.join(f'{x:.6E}' for x in he[:2]) + '  ! zone 99\n   \n   ' + '   '.join(f'{x:.6E}' for x in he[2:])
        p = tmp_path / 'SN_HYDRO_DATA'
        p.write_text('\n'.join(lines) + '\n')
        hydro = pv.parse_sn_hydro(p)
        np.testing.assert_allclose(hydro.mass_fractions['HE'], he[::-1], rtol=1e-5)
        np.testing.assert_allclose(hydro.mass_fractions['HYD'], (1.0 - he)[::-1], rtol=1e-5)

def _hydro_dataclass(nd=6, ed_min=1.0):
    r = np.linspace(2.4e4, 1.0e4, nd)
    return pp.HydroData(
        nd=nd, time_days=1.0, lum_lsun=1e6, r_inner_cm=1e14, t_inner_k=2e4,
        radius=r, velocity=np.linspace(6000.0, 1000.0, nd),
        sigma=np.zeros(nd), temperature=np.linspace(0.6, 2.0, nd),
        density=np.linspace(1e-12, 1e-9, nd),
        ed=np.geomspace(ed_min, 1e12, nd),
        natom=np.linspace(1e10, 1e14, nd), kappa=np.full(nd, 0.34),
        mass_fracs={'HYD': np.linspace(0.9, 0.4, nd),
                    'HE': np.linspace(0.05, 0.55, nd),
                    'OXY': np.linspace(0.05, 0.05, nd)},
        raw_path='')

def _synth_atomic(name='HI', n_levels=3, ionization_en=109678.764):
    stat_wt = np.array([2.0 * (i + 1) ** 2 for i in range(n_levels)])
    energy = np.zeros(n_levels)
    for i in range(1, n_levels):
        energy[i] = ionization_en * (1.0 - 1.0 / (i + 1) ** 2)
    fedge = (ionization_en - energy) * pp.SPEED_LIGHT * 1.0e-15
    return pp.AtomicData(name=name, n_levels=n_levels,
                         ionization_energy=ionization_en, stat_wt=stat_wt,
                         energy=energy, fedge=fedge, transitions=[])

def _write_foscdat(path, n_levels=3, ionization_en=109678.764):
    lines = [
        '12-OCT-1995                                      !Date',
        f'{n_levels}                                               !Number of energy levels',
        f'{ionization_en:.4f}                                      !Ionization energy',
        '1.0                                              !Screened nuclear charge',
        '0                                                !Number of transitions',
        '',
    ]
    for i in range(n_levels):
        g = 2.0 * (i + 1) ** 2
        e_cm = ionization_en * (1.0 - 1.0 / (i + 1) ** 2) if i > 0 else 0.0
        freq = (ionization_en - e_cm) * pp.SPEED_LIGHT * 1.0e-15
        lines.append(f'{i+1:3d}        {g:10.4f}      {e_cm:14.8f}     {freq:10.6f}     {i+1}')
    Path(path).write_text('\n'.join(lines) + '\n')

def _make_model_dir(root, vel_law='11', nd=6, ed_min=1.0, gammas_ed=1e20):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    ed = np.geomspace(ed_min, 1e12, nd)
    lines = _hydro_lines(nd=nd, ed=ed)
    (root / 'SN_HYDRO_DATA').write_text('\n'.join(lines) + '\n')
    (root / 'MODEL_SPEC').write_text('\n'.join([
        f'{nd} [ND]', '8 [NC]', '14 [NP]', '1 [NUM_BNDS]', '!',
        '3,3,3 [HI_ISF]', '3,3,3 [HeI_ISF]',
    ]) + '\n')
    (root / 'VADAT').write_text(
        f'{vel_law} [VEL_LAW]\nT [LTE_EST]\nF [LIN_INT]\nTR [DC_METH]\n'
        '-9.9D-01 [HYD/X]\n-9.9D-01 [HE/X]\n-9.9D-01 [OXY/X]\n')
    for name in ('HI', 'HeI'):
        _write_foscdat(root / f'{name}_F_OSCDAT')
    gam = ['!Generated test', f'!Number of depth points\n{nd}',
           '!Electron density',
           '   ' + '   '.join(f'{gammas_ed:10.4E}' for _ in range(nd)),
           '!Radius (10^10cm)',
           '   ' + '   '.join(f'{x:10.4E}' for x in np.linspace(2.4e4, 1.0e4, nd))]
    (root / 'GAMMAS_IN').write_text('\n'.join(gam) + '\n')
    return root

class TestHydroReferenceAbundances:

    def test_vel_law_11_uses_depth_dependent_fracs(self):
        hydro = _hydro_dataclass()
        hydro.mass_fracs['HE'][2] = 0.0
        inputs = pp.ModelInputs(
            hydro=hydro, species=[], vadat_abundances={'HYD': 0.01, 'HE': 0.99},
            vadat_flags={'VEL_LAW': '11'}, atom_data={})
        conv = pp.CmfgenConventions(strategy='LTE')
        ab = pp.build_abundances(inputs, conv)
        assert set(ab) == set(hydro.mass_fracs)
        for e, fracs in hydro.mass_fracs.items():
            np.testing.assert_array_equal(ab[e], fracs)
            assert ab[e] is not fracs
        assert ab['HE'][2] == 0.0

    def test_nonhydro_legacy_vadat_preferred(self):
        hydro = _hydro_dataclass()
        inputs = pp.ModelInputs(
            hydro=hydro, species=[], vadat_abundances={'HYD': 0.01},
            vadat_flags={'VEL_LAW': '12'}, atom_data={})
        conv = pp.CmfgenConventions(strategy='LTE')
        ab = pp.build_abundances(inputs, conv)
        np.testing.assert_allclose(ab['HYD'], np.full(hydro.nd, 0.01))
        np.testing.assert_array_equal(ab['HE'], hydro.mass_fracs['HE'])

    def test_vel_law_11_without_fracs_raises(self):
        hydro = _hydro_dataclass()
        hydro.mass_fracs = {}
        inputs = pp.ModelInputs(
            hydro=hydro, species=[], vadat_abundances={'HYD': 0.5},
            vadat_flags={'VEL_LAW': '11'}, atom_data={})
        with pytest.raises(ValueError):
            pp.build_abundances(inputs, pp.CmfgenConventions())

class TestGenerateHydroReference:

    def test_source_ed_preserved_and_source_untouched(self, tmp_path, monkeypatch):
        mdir = _make_model_dir(tmp_path / 'model', vel_law='11', ed_min=1.0)
        out = tmp_path / 'proxy_out'
        inputs = pp.load_model(str(mdir))
        src_ed = inputs.hydro.ed.copy()
        assert src_ed.min() == pytest.approx(1.0)

        captured = []
        orig = pp.compute_all_element_ionizations
        def spy(*a, **kw):
            captured.append(kw.get('floor_ne'))
            return orig(*a, **kw)
        monkeypatch.setattr(pp, 'compute_all_element_ionizations', spy)

        snapshots = {f.name: f.read_bytes() for f in mdir.iterdir() if f.is_file()}
        conv = pp.CmfgenConventions(strategy='LTE', ne_max_iter=1,
                                    floor_ne=10.0, smooth=False)
        result = pp.generate(inputs, conv, str(out))

        assert captured and all(f == pytest.approx(1.0) for f in captured)

        in_files = sorted(out.glob('*_IN'))
        assert in_files, 'no _IN files written'
        for f in in_files:
            data = pp.read_in_file(str(f))
            assert data is not None
            np.testing.assert_allclose(data.ed, src_ed, rtol=1e-8)

        for f in mdir.iterdir():
            if f.is_file():
                assert f.read_bytes() == snapshots[f.name], f'{f.name} modified'

        ne_log = result['strategy_log']['ne_iteration']
        assert ne_log['enabled'] is False
        assert ne_log['converged'] is None
        assert ne_log['reference'] == 'SN_HYDRO_DATA'
        assert 'disabled' in ne_log['note']
        assert result['strategy_log']['tau_trans'] == conv.tau_trans

    def test_iterating_ne_on_hydro_reference_raises_before_makedirs(self, tmp_path):
        mdir = _make_model_dir(tmp_path / 'model', vel_law='11')
        out = tmp_path / 'proxy_out_ne2'
        inputs = pp.load_model(str(mdir))
        conv = pp.CmfgenConventions(strategy='LTE', ne_max_iter=5, smooth=False)
        with pytest.raises(ValueError, match='VEL_LAW=11'):
            pp.generate(inputs, conv, str(out))
        assert not out.exists()

class TestInputCompletenessGate:

    def test_missing_atomic_data_fails_before_output(self, tmp_path):
        mdir = _make_model_dir(tmp_path / 'model', vel_law='12')
        (mdir / 'HI_F_OSCDAT').unlink()
        out = tmp_path / 'proxy_out_missing'
        inputs = pp.load_model(str(mdir))
        conv = pp.CmfgenConventions(strategy='LTE', ne_max_iter=1, smooth=False)
        with pytest.raises(ValueError, match='Missing atomic data'):
            pp.generate(inputs, conv, str(out))
        assert not out.exists()

    def test_no_active_species_fails_before_output(self, tmp_path):
        mdir = _make_model_dir(tmp_path / 'model', vel_law='12')
        (mdir / 'MODEL_SPEC').write_text('\n'.join([
            '6 [ND]', '8 [NC]', '14 [NP]', '1 [NUM_BNDS]', '!',
            '3,3,0 [HI_ISF]', '3,3,0 [HeI_ISF]',
        ]) + '\n')
        out = tmp_path / 'proxy_out_noact'
        inputs = pp.load_model(str(mdir))
        conv = pp.CmfgenConventions(strategy='LTE', ne_max_iter=1, smooth=False)
        with pytest.raises(ValueError, match='No active species'):
            pp.generate(inputs, conv, str(out))
        assert not out.exists()

    def test_failed_validation_blocks_vadat_patch(self, tmp_path, monkeypatch):
        mdir = _make_model_dir(tmp_path / 'model', vel_law='12')
        out = tmp_path / 'proxy_out_val'
        out.mkdir()
        vadat_copy = out / 'VADAT'
        vadat_copy.write_bytes((mdir / 'VADAT').read_bytes())
        before = vadat_copy.read_bytes()
        monkeypatch.setattr(pp, 'validate_run',
                            lambda *a, **kw: {'issues': ['forced'], 'n_issues': 1,
                                              'ne_ratio_max': 0.0})
        inputs = pp.load_model(str(mdir))
        conv = pp.CmfgenConventions(strategy='LTE', ne_max_iter=1, smooth=False)
        pp.generate(inputs, conv, str(out))
        assert vadat_copy.read_bytes() == before
        rc = pp.main([str(mdir), '--output-dir', str(tmp_path / 'o2'),
                      '--strategy', 'LTE', '--ne-max-iter', '1'])
        assert rc == 1


class TestSelfConsistentED:

    def test_written_ed_matches_last_ionization_input(self, tmp_path, monkeypatch):
        mdir = _make_model_dir(tmp_path / 'model', vel_law='12')
        out = tmp_path / 'proxy_out'
        inputs = pp.load_model(str(mdir))

        seen_ed = []
        orig = pp.compute_all_element_ionizations
        def spy(*a, **kw):
            seen_ed.append(inputs.hydro.ed.copy())
            return orig(*a, **kw)
        monkeypatch.setattr(pp, 'compute_all_element_ionizations', spy)

        conv = pp.CmfgenConventions(strategy='LTE', ne_max_iter=5,
                                    ne_tol=1e-12, smooth=False)
        pp.generate(inputs, conv, str(out))

        assert len(seen_ed) >= 1
        last_ed = seen_ed[-1]
        for f in out.glob('*_IN'):
            data = pp.read_in_file(str(f))
            np.testing.assert_allclose(data.ed, last_ed, rtol=1e-8)

def _writer_kwargs(r_order, src_order):
    n = 4
    r = np.array([4.0, 3.0, 2.0, 1.0]) * 1e14
    if r_order == 'asc':
        r = r[::-1]
    src_r = np.array([4.0, 3.0, 2.0, 1.0]) * 1e14
    src_ne = np.array([40.0, 30.0, 20.0, 10.0])
    if src_order == 'asc':
        src_r, src_ne = src_r[::-1], src_ne[::-1]
    return dict(
        r_cm=r.copy(), v_kms=np.full(n, 1000.0), rho=np.full(n, 1e-10),
        T_K=np.full(n, 4000.0), sigma=np.zeros(n), kappa=np.full(n, 0.34),
        time_days=30.0, mass_fractions={'H': np.ones(n)}, lum_lsun=1e6,
        t_floor=0.0, inc_rad_decays=False,
        ne_stella=src_ne, ne_stella_r=src_r)

class TestWriterElectronDensity:

    @pytest.mark.parametrize('r_order', ['asc', 'desc'])
    @pytest.mark.parametrize('src_order', ['asc', 'desc'])
    def test_source_ordering_independent(self, tmp_path, r_order, src_order):
        out = tmp_path / f'SN_HYDRO_DATA_{r_order}_{src_order}'
        kw = _writer_kwargs(r_order, src_order)
        src_r = kw['ne_stella_r'].copy()
        src_ne = kw['ne_stella'].copy()
        cvt.write_sn_hydro_data(str(out), **kw)
        d = read_sn_hydro_data(out)
        order = np.argsort(src_r)
        expected = np.interp(d.r_1e10cm * 1e10, src_r[order], src_ne[order])
        np.testing.assert_allclose(d.electron_density, expected, rtol=1e-6)

    def test_unpaired_source_args_fail(self, tmp_path):
        kw = _writer_kwargs('desc', 'desc')
        kw['ne_stella_r'] = None
        with pytest.raises(ValueError):
            cvt.write_sn_hydro_data(str(tmp_path / 'o7'), **kw)
        kw = _writer_kwargs('desc', 'desc')
        kw['ne_stella'] = None
        with pytest.raises(ValueError):
            cvt.write_sn_hydro_data(str(tmp_path / 'o8'), **kw)

    def test_mismatched_lengths_fail(self, tmp_path):
        kw = _writer_kwargs('desc', 'desc')
        kw['ne_stella'] = kw['ne_stella'][:-1]
        with pytest.raises(ValueError):
            cvt.write_sn_hydro_data(str(tmp_path / 'o1'), **kw)

    def test_duplicate_source_radius_fails(self, tmp_path):
        kw = _writer_kwargs('desc', 'desc')
        kw['ne_stella_r'] = kw['ne_stella_r'].copy()
        kw['ne_stella_r'][1] = kw['ne_stella_r'][0]
        with pytest.raises(ValueError):
            cvt.write_sn_hydro_data(str(tmp_path / 'o2'), **kw)

    def test_out_of_range_target_fails(self, tmp_path):
        kw = _writer_kwargs('desc', 'desc')
        kw['r_cm'] = np.array([5.0, 3.0, 2.0, 1.0]) * 1e14
        with pytest.raises(ValueError):
            cvt.write_sn_hydro_data(str(tmp_path / 'o3'), **kw)

    def test_writer_rejects_bad_inputs(self, tmp_path):
        kw = _writer_kwargs('desc', 'desc')
        kw['time_days'] = -5.0
        with pytest.raises(ValueError):
            cvt.write_sn_hydro_data(str(tmp_path / 'o4'), **kw)
        kw = _writer_kwargs('desc', 'desc')
        kw['rho'] = np.array([1e-10, 0.0, 1e-10, 1e-10])
        with pytest.raises(ValueError):
            cvt.write_sn_hydro_data(str(tmp_path / 'o5'), **kw)
        kw = _writer_kwargs('desc', 'desc')
        kw['mass_fractions'] = {'H': np.array([1.0, 1.1, 0.9, 1.0])}
        with pytest.raises(ValueError):
            cvt.write_sn_hydro_data(str(tmp_path / 'o6'), **kw)

    def test_boundary_metadata_matches_written_file(self, tmp_path):
        out = tmp_path / 'SN_HYDRO_DATA'
        n = 4
        cvt.write_sn_hydro_data(
            str(out), r_cm=np.array([4.0, 3.0, 2.0, 1.0]) * 1e14,
            v_kms=np.full(n, 1000.0), rho=np.full(n, 1e-10),
            T_K=np.array([2000.0, 2000.0, 4000.0, 5000.0]),
            sigma=np.zeros(n), kappa=np.full(n, 0.34),
            time_days=30.0, mass_fractions={'H': np.ones(n)}, lum_lsun=1e6,
            t_floor=3000.0, inc_rad_decays=False,
            r_inner_cm=0.9e14, T_inner_K=2000.0, tau_inner=0.5, tau_inner_chi=0.4)
        d = read_sn_hydro_data(out)
        hdr_r = float(d.headers['Inner boundary radius [cm]'])
        hdr_t = float(d.headers['Inner boundary temperature [K]'])
        hdr_tau = float(d.headers['Inner boundary Rosseland tau'])
        hdr_tau_chi = float(d.headers['Inner boundary chi-column tau'])
        assert hdr_r == pytest.approx(d.r_1e10cm[-1] * 1e10, rel=1e-6)
        assert hdr_t == pytest.approx(d.T_1e4K[-1] * 1e4, rel=1e-6)
        tau_file = cvt.compute_kappa_tau(d.r_1e10cm * 1e10, d.density_gcm3, d.kappa)
        assert hdr_tau == pytest.approx(tau_file[-1], rel=1e-6)
        assert hdr_tau_chi == pytest.approx(tau_file[-1], rel=1e-6)
        assert float(d.headers['Source reference boundary radius [cm]']) == pytest.approx(0.9e14)
        assert float(d.headers['Source reference boundary temperature [K]']) == pytest.approx(2000.0)
        assert float(d.headers['Pre-write reference boundary tau']) == pytest.approx(0.5)
        assert float(d.headers['Pre-write chi-column tau']) == pytest.approx(0.4)

class TestSubmitPlumbing:

    def _stub_base(self, tmp_path):
        base = tmp_path / 'stub_base'
        base.mkdir()
        args_file = tmp_path / 'stub_args.txt'
        (base / 'submit_cmfgen.sh').write_text(
            '#!/bin/bash\n'
            f'printf "%s\\n" "$@" > "{args_file}"\n'
            'printf "CMFGEN Job ID: 111\\nOBS Job ID: 222\\n" > submit.log\n')
        (base / 'submit_cmfgen.sh').chmod(0o755)
        return base, args_file

    def _run_jobs(self, tmp_path, extra):
        base, args_file = self._stub_base(tmp_path)
        cmfgen_dir = tmp_path / 'cmfgen_dir'
        cmfgen_dir.mkdir(exist_ok=True)
        env = dict(os.environ, CMFGEN_BASE_DIR=str(base))
        res = subprocess.run(
            ['bash', str(SCRIPTS / 'submit_cmfgen_jobs.sh'),
             'M', str(tmp_path / 'model'), str(cmfgen_dir), '0'] + extra,
            capture_output=True, text=True, env=env)
        args = args_file.read_text() if args_file.exists() else ''
        return res, args

    def test_non_local_reversals_flag_forwarded(self, tmp_path):
        res, args = self._run_jobs(tmp_path, ['--proxy-no-non-local-reversals=1'])
        assert res.returncode == 0, res.stderr
        assert '--proxy-no-non-local-reversals=1' in args

    def test_non_local_reversals_zero_dropped(self, tmp_path):
        res, args = self._run_jobs(tmp_path, ['--proxy-no-non-local-reversals=0'])
        assert res.returncode == 0, res.stderr
        assert '--proxy-no-non-local-reversals=1' not in args

    def test_non_local_reversals_omitted(self, tmp_path):
        res, args = self._run_jobs(tmp_path, [])
        assert res.returncode == 0, res.stderr
        assert 'no-non-local-reversals' not in args

    def test_submit_cmfgen_bool_parser_validates(self):
        src = (REPO / 'cmfgen_base' / 'submit_cmfgen.sh').read_text()
        m = re.search(r'--proxy-no-non-local-reversals=\*\)(.*?shift ;;)', src, re.S)
        assert m is not None
        block = m.group(1)
        assert '0) PROXY_NO_NON_LOCAL_REVERSALS=""' in block
        assert '1) PROXY_NO_NON_LOCAL_REVERSALS="1"' in block
        assert 'exit 2' in block

    def test_recognised_proxy_params_forwarded(self, tmp_path):
        res, args = self._run_jobs(tmp_path, [
            '--proxy-ne-max-iter=1', '--proxy-tau-trans=1.5',
            '--proxy-ionization-mode=lte', '--proxy-v2-strategy=GREY'])
        assert res.returncode == 0, res.stderr
        for frag in ('--proxy-ne-max-iter=1', '--proxy-tau-trans=1.5',
                     '--proxy-ionization-mode=lte', '--proxy-v2-strategy=GREY'):
            assert frag in args

    def test_unknown_proxy_flag_rejected(self, tmp_path):
        res, args = self._run_jobs(tmp_path, ['--proxy-typo=x'])
        assert res.returncode == 2
        assert 'unknown proxy option' in res.stderr
        assert args == ''

    def test_v2_proxy_command_has_no_or_true(self):
        src = (REPO / 'cmfgen_base' / 'submit_cmfgen.sh').read_text().splitlines()
        v2 = [ln for ln in src if 'proxy_populations_v2.py' in ln and '>>' in ln]
        assert len(v2) == 1
        assert '|| true' not in v2[0]

    def test_proxy_failure_terminates_block(self, tmp_path):
        stub = tmp_path / 'fake_proxy.py'
        stub.write_text('import sys; sys.exit(7)\n')
        block = tmp_path / 'block.sh'
        marker = tmp_path / 'launched'
        block.write_text(
            'set -euo pipefail\n'
            f'python "{stub}" >> "batch.log"\n'
            f'echo CMFGEN_LAUNCHED > "{marker}"\n')
        res = subprocess.run(['bash', str(block)], cwd=tmp_path,
                             capture_output=True, text=True)
        assert res.returncode == 7
        assert not marker.exists()


import run_model_grid as rmg


class TestInlineConversionGate:
    """run_inline_cmfgen_conversion must fail on nonzero converter exit."""

    def _setup(self, tmp_path, monkeypatch, returncode, stdout, stderr):
        model_dir = tmp_path / 'model'
        cmfgen_dir = model_dir / 'cmfgen' / '5.00days'
        cmfgen_dir.mkdir(parents=True)
        sn_hydro = cmfgen_dir / 'SN_HYDRO_DATA'
        sn_hydro.write_bytes(b'OLD')
        config = {'output': {'model_dir': str(model_dir)}}
        calls = []

        def fake_run(cmd, **kw):
            calls.append(cmd)
            sn_hydro.write_bytes(b'NEW')
            return subprocess.CompletedProcess(cmd, returncode,
                                               stdout=stdout, stderr=stderr)

        monkeypatch.setattr(rmg, 'build_cmfgen_conversion_cmd',
                            lambda *a, **k: 'fake converter')
        monkeypatch.setattr(rmg.subprocess, 'run', fake_run)
        monkeypatch.setattr(rmg, 'generate_grey_scl_fac_in', lambda p: None)
        monkeypatch.setattr(rmg, 'check_sn_hydro_sigma', lambda *a, **k: None)
        return sn_hydro, config, calls

    def test_nonzero_exit_fails_and_preserves_artifacts(self, tmp_path,
                                                      monkeypatch):
        sn_hydro, config, calls = self._setup(
            tmp_path, monkeypatch, 2, 'Post-write audit FAIL', 'failure')
        ok = rmg.run_inline_cmfgen_conversion(config, tmp_path, tmp_path, 5.0)
        assert ok is False
        assert sn_hydro.read_bytes() == b'NEW'
        attempts = list(sn_hydro.parent.glob('conversion_attempt_*'))
        assert len(attempts) == 1
        attempt = attempts[0]
        assert (attempt / 'SN_HYDRO_DATA.previous').read_bytes() == b'OLD'
        assert (attempt / 'command.txt').read_text() == 'fake converter\n'
        assert (attempt / 'stdout.log').read_text() == 'Post-write audit FAIL'
        assert (attempt / 'stderr.log').read_text() == 'failure'
        assert calls == ['fake converter']

    def test_success_saves_stdout(self, tmp_path, monkeypatch):
        sn_hydro, config, calls = self._setup(
            tmp_path, monkeypatch, 0, 'audit ok', '')
        ok = rmg.run_inline_cmfgen_conversion(config, tmp_path, tmp_path, 5.0)
        assert ok is True
        attempts = list(sn_hydro.parent.glob('conversion_attempt_*'))
        assert len(attempts) == 1
        assert (attempts[0] / 'stdout.log').read_text() == 'audit ok'
        assert (attempts[0] / 'SN_HYDRO_DATA.previous').read_bytes() == b'OLD'
        assert calls == ['fake converter']


class TestPopSpeciesFloor:

    def test_proxy_di_floored_where_species_absent(self, tmp_path):
        nd = 6
        mdir = _make_model_dir(tmp_path / 'model', vel_law='11', nd=nd)
        fracs = {'HYD': np.array([0.0, 0.0, 0.3, 0.5, 0.6, 0.7]),
                 'HE': np.array([0.3, 0.3, 0.3, 0.3, 0.0, 0.0])}
        _write_hydro(mdir / 'SN_HYDRO_DATA', nd=nd, fracs=fracs)
        out = tmp_path / 'proxy_out_floor'
        inputs = pp.load_model(str(mdir))
        conv = pp.CmfgenConventions(strategy='LTE', ne_max_iter=1, smooth=False)
        result = pp.generate(inputs, conv, str(out))
        in_files = sorted(out.glob('*_IN'))
        assert in_files
        for f in in_files:
            data = pp.read_in_file(str(f))
            assert np.all(np.isfinite(data.di))
            assert np.all(data.di > 0), f'{f.name} has DI<=0'
            assert np.all(data.di >= 1e-300)
        floored = result['strategy_log']['pop_species_floor']
        assert floored['value_cm3'] == 1e-20
        assert floored['species_depths_floored']['HI'] == 2
        assert floored['species_depths_floored']['HeI'] == 2

    def test_write_in_file_rejects_nonpositive_di(self, tmp_path):
        nd = 4
        hydro = _hydro_dataclass(nd=nd)
        nlev = 2
        di = np.full(nd, 1e10)
        di[1] = 0.0
        pop = pp.SpeciesPopulations(
            name='HI', elem_label='HYD', ion_stage=1, nlev=nlev,
            di=di, dc=np.ones((nlev, nd)),
            t_write=hydro.temperature.copy(),
            target_pop=np.ones((nlev, nd)), lte_pop=np.ones((nlev, nd)),
            ion_frac_current=np.ones(nd), ion_frac_next=np.ones(nd),
            omitted_stage_frac=np.ones(nd),
            floored_mask=np.zeros((nlev, nd), dtype=bool),
            capped_mask=np.zeros((nlev, nd), dtype=bool),
            abundance_residual=np.zeros(nd), smooth_diag={})
        target = tmp_path / 'HI_IN'
        with pytest.raises(ValueError, match='DI must be finite and > 0'):
            pp.write_in_file(str(target), pop, hydro,
                             pp.CmfgenConventions())
        assert not target.exists()
        assert not Path(str(target) + '.tmp').exists()

    def test_cleanup_keeps_conversion_attempt_dirs(self, tmp_path):
        src = (REPO / 'cmfgen_base' / 'submit_cmfgen.sh').read_text().splitlines()
        start = next(i for i, ln in enumerate(src)
                     if ln.strip() == 'for item in *; do')
        indent = len(src[start]) - len(src[start].lstrip())
        end = next(i for i in range(start, len(src))
                   if src[i].strip() == 'done'
                   and len(src[i]) - len(src[i].lstrip()) == indent)
        loop = '\n'.join(src[start:end + 1])
        assert 'conversion_attempt_*' in loop
        (tmp_path / 'foo').write_text('x')
        (tmp_path / 'conversion_attempt_abc').mkdir()
        (tmp_path / 'conversion_attempt_abc' / 'x').write_text('y')
        (tmp_path / 'VADAT').write_text('v')
        script = (
            'set -e\n'
            'FILES_TO_KEEP="VADAT SN_HYDRO_DATA"\n'
            f'cd "{tmp_path}"\n' + loop + '\n')
        res = subprocess.run(['bash', '-c', script], capture_output=True,
                             text=True)
        assert res.returncode == 0, res.stderr
        assert not (tmp_path / 'foo').exists()
        assert (tmp_path / 'conversion_attempt_abc' / 'x').exists()
        assert (tmp_path / 'VADAT').exists()
