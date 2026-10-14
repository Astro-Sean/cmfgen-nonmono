#!/usr/bin/env python3
"""Mathematical verification of the JFNK+FIX_T_AUTO invariants.

Drives bin/se_driver (which links the production jfnk_solver_v1
newton_hook_v1 / gmres_m_v1 and the surrogate replica of the CMFGEN
interface) on the generated systems.  Evidence for the audit report.
"""

import os
import subprocess
import sys
import tempfile

import numpy as np

DRV = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "bin", "se_driver")
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "systems")
NT, ND = 37, 4
TROW = NT - 1          # 1-based T row inside a depth block
TIDX = [(k - 1) * NT + TROW - 1 for k in range(1, ND + 1)]  # 0-based


def run(*args):
    r = subprocess.run([DRV, *map(str, args)], capture_output=True,
                       text=True)
    return r.stdout


def vec_file(v):
    fd, p = tempfile.mkstemp(suffix=".txt")
    with os.fdopen(fd, "w") as f:
        f.write("\n".join("%.17e" % x for x in v))
    return p


def load_vec(p):
    return np.loadtxt(p)


def sejvp(spec, x, v, eps, tauto=1):
    vf, xf = vec_file(v), vec_file(x)
    out = run("sejvp", os.path.join(ROOT, spec), vf, eps, xf, tauto)
    os.unlink(vf)
    os.unlink(xf)
    d = {}
    jv_fd, jv_ex = [], []
    for ln in out.splitlines():
        if ln.startswith("JVFDNORM"):
            d["fd"] = float(ln.split()[1])
        elif ln.startswith("JVEXNORM"):
            d["ex"] = float(ln.split()[1])
        elif ln.startswith("RELDIFF"):
            d["rel"] = float(ln.split()[1])
        elif ln.startswith("JV "):
            p = ln.split()
            jv_fd.append(float(p[2]))
            jv_ex.append(float(p[3]))
    return d, np.array(jv_fd), np.array(jv_ex)


def sef(spec, x, tauto=1):
    xf = vec_file(x)
    out = run("sef", os.path.join(ROOT, spec), xf, tauto)
    os.unlink(xf)
    f = []
    for ln in out.splitlines():
        p = ln.split()
        if len(p) == 3 and p[0] == "F":
            f.append(float(p[2]))
    return np.array(f)


def main():
    rng = np.random.default_rng(20261013)
    n = NT * ND
    xstar = np.loadtxt(os.path.join(ROOT, "shock", "xstar.txt"))
    states = {
        "x=0 (STELLA-like, |RE| large)": np.zeros(n),
        "x=x*/2 (mid)": 0.5 * xstar,
        "x=x* (converged)": xstar,
    }

    print("=" * 72)
    print("INVARIANT 1: FD-JVP vs analytic J*v (shock spec, tauto=1)")
    print("=" * 72)
    epss = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
    for sname, x in states.items():
        print(f"\n  state {sname}")
        rel_errs = {e: [] for e in epss}
        for t in range(5):
            v = rng.standard_normal(n)
            for e in epss:
                d, _, _ = sejvp("shock", x, v, e, tauto=1)
                rel_errs[e].append(d["rel"])
        for e in epss:
            errs = rel_errs[e]
            print(f"    eps={e:8.0e}  reldiff min/med/max = "
                  f"{min(errs):.3e} {np.median(errs):.3e} {max(errs):.3e}")

    print("\n  -- sub-check 1a: v_T only (pop comps = 0)")
    x = np.zeros(n)
    vt = np.zeros(n)
    vt[TIDX] = 1.0
    d, jfd, jex = sejvp("shock", x, vt, 1e-6, tauto=1)
    pop_ex = np.delete(jex, TIDX)
    print(f"    |Jv_ex|={d['ex']:.4e}  max|pop rows|={np.abs(pop_ex).max():.4e}"
          f"  T rows={np.abs(jex[TIDX]).max():.4e}  reldiff={d['rel']:.3e}")
    print("    -> dR_SE/dT nonzero:", np.abs(pop_ex).max() > 0,
          " dR_E/dT nonzero:", np.abs(jex[TIDX]).max() > 0)

    print("\n  -- sub-check 1b: v_n only (T comps = 0)")
    vn = rng.standard_normal(n)
    vn[TIDX] = 0.0
    d, jfd, jex = sejvp("shock", x, vn, 1e-6, tauto=1)
    print(f"    |Jv_ex|={d['ex']:.4e}  max|T rows|={np.abs(jex[TIDX]).max():.4e}"
          f"  reldiff={d['rel']:.3e}")
    print("    -> dR_E/dn nonzero:", np.abs(jex[TIDX]).max() > 0)

    print("\n" + "=" * 72)
    print("INVARIANT 2/7: damped admission, stored vs fresh F, trajectory")
    print("=" * 72)
    out = run("setauto", os.path.join(ROOT, "shock"),
              "DDDDDDDDDDDD", "0.05", "0.1", "0.1", "0.5", "2")
    print("   adm   INFO   eps_T      rung alpha    XT          "
          "|F|stored vs fresh (reldiff)")
    xt_hist, eps_hist = [], []
    for ln in out.splitlines():
        p = ln.split()
        if ln.startswith("AUTO "):
            eps_hist.append(float(p[9]))
            xt_hist.append(None)
        elif ln.startswith("FCHK "):
            i = p[1]
            st, fr = float(p[3]), float(p[5])
            rel = abs(st - fr) / max(fr, 1e-300)
            print(f"   {i:>3}  stored={st:.10e} fresh={fr:.10e} "
                  f"reldiff={rel:.3e}")
        elif ln.startswith("XTALL "):
            xt_hist[int(p[1]) - 1] = np.array(
                [float(t) for t in p[2:]])
    # continuity: per-depth |xT_k(d) - xT_{k-1}(d)| <= eps_T,k and at
    # least one depth attains it; jumps must equal the committed step.
    print("\n   trajectory continuity (per-depth):")
    ok = True
    xprev = np.zeros(ND)
    for k, xt in enumerate(xt_hist):
        if xt is None:
            continue
        d = np.abs(xt - xprev)
        # eps_T prints at 4 sig figs (ES10.3): compare at print precision
        good = abs(d.max() - eps_hist[k]) < 6e-4 * max(1.0, eps_hist[k])
        ok &= good
        ds = np.array2string(d, precision=3, separator=",")
        print(f"   adm {k+1:2d}: per-depth dxT={ds}  max={d.max():.4e} "
              f"eps_T={eps_hist[k]:.4e}  {'OK' if good else 'MISMATCH'}")
        xprev = xt
    print("   continuous:", ok)

    print("\n" + "=" * 72)
    print("INVARIANT 5: R_E reads the current T (no stale state)")
    print("=" * 72)
    f0 = sef("shock", np.zeros(n), tauto=1)
    f0b = sef("shock", np.zeros(n), tauto=1)
    print(f"   determinism: max|F-F|={np.abs(f0 - f0b).max():.3e}")
    fstar = sef("shock", xstar, tauto=1)
    print(f"   at x*: |F|={np.linalg.norm(fstar):.4e}  "
          f"|RE|(T rows)={np.linalg.norm(fstar[TIDX]):.4e}")
    xp = np.zeros(n)
    xp[TIDX[0]] = 0.01          # +1% in log-T at depth 1
    fp = sef("shock", xp, tauto=1)
    ch = fp - f0
    ch_t = np.abs(ch[TIDX])
    ch_pop = np.delete(np.abs(ch), TIDX)
    print(f"   +1% log-T at depth1: max|dF|={np.abs(ch).max():.4e}")
    print(f"     T-row responses: {ch_t}")
    print(f"     pop-row response at depth1 rows 0..{NT-3}: "
          f"{np.abs(ch[:NT-2]).max():.4e}")
    print("     R_E responds to T perturbation:", ch_t[0] > 0)

    print("\n" + "=" * 72)
    print("INVARIANT 3/4: GMRES on shock, T rows live, PC scale on/off")
    print("=" * 72)
    for mode, lbl in [(0, "T masked (baseline)"),
                      (1, "T live, PC unscaled"),
                      (2, "T live, PC scaled a=0.05")]:
        out = run("segmres", os.path.join(ROOT, "shock"),
                  "20", "30", "6", "1e-2", "1", "0", "0",
                  "1e-2", "1", str(mode), "0.05")
        last = [ln for ln in out.splitlines()
                if ln.strip() and ln.split()[0].isdigit()]
        if last:
            k, its, res, xinf, info = last[-1].split()
            first = last[0].split()
            print(f"   {lbl:28s}: res {float(first[2]):.3e} -> "
                  f"{float(res):.3e} in {its} its (info={info})")

    print("\n" + "=" * 72)
    print("INVARIANT 6: damp index map (k-1)*NT+j  vs  expected", TIDX)
    print("=" * 72)
    out = run("sedamp", os.path.join(ROOT, "shock"), "0.05", "0.5",
              "1", "4")
    hit = []
    for ln in out.splitlines():
        p = ln.split()
        if ln.startswith("DAMPED"):
            hit.append(int(p[1]) - 1)
            print("   ", ln.strip())
    print("   matches expected T indices:", sorted(hit) == sorted(TIDX))
    out = run("sedamp", os.path.join(ROOT, "shock"), "0.05", "0.5",
              "2", "3")
    hit2 = [int(ln.split()[1]) - 1 for ln in out.splitlines()
            if ln.startswith("DAMPED")]
    print("   partial domain [2,3] hits:", hit2,
          "expected:", [TIDX[1], TIDX[2]])


if __name__ == "__main__":
    main()
