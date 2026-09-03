# Audit: STELLA Radiation Field → CMFGEN EDDFACTOR Initialisation

## 1. STELLA .swd Current Output

### Source File
`/ptmp/sbrennan/stella/OpenStella_conv/strad/lbalsw.f`, subroutine `LOSSEN`, lines 1038-1055.

### Variables Currently Written per Zone (unit 81)

| Column | Variable | Physical Definition | Units | Frame |
|--------|----------|---------------------|-------|-------|
| 1 | `t_day` | Time since explosion | days | — |
| 2 | `km` | Zone index (1=innermost, NZON=outermost) | — | — |
| 3 | `log10(AMPR)` | Enclosed mass (external layers) | Msun | — |
| 4 | `log10(UR*Ry(Km))` | Radius | cm | — |
| 5 | `Uy*1e6/(UTIME*CRAP)` | Velocity | km/s | lab |
| 6 | `log10(UTP*Ty)` | Gas temperature | K | comoving |
| 7 | `log10(UTP*TpRAD)` | Radiation temperature | K | comoving |
| 8 | `PLLOG` | log10 density | g/cm^3 (internal) | — |
| 9 | `PLOG` | log10 pressure | cgs (internal) | — |
| 10 | `QVLOG` | log10 pseudo-viscosity | internal | — |
| 11 | `log10(eng)` | Radiation energy density | internal | comoving |
| 12 | `Flum*1e-40` | Luminosity | 1e40 erg/s | lab |
| 13 | `WRKX(Km)` | Rosseland opacity | internal | — |
| 14 | `barionn` | Baryon number density | cm^-3 | — |
| 15 | `Xe*barionn` | Electron number density | cm^-3 | — |

Composition (15 elements) is written after these columns.

### Key Limitation
The .swd contains **only bolometric** radiation quantities (T_rad, eng, Flum).
**No frequency-dependent information is present.**

---

## 2. STELLA Internal Radiation Variables (Available but Not Output)

### Frequency-Group Structure
STELLA uses `NFREQ = 100` frequency groups (parameter `NFREQ=96+4`).

### Frequency Grid (COMMON/BAND/)
| Variable | Dimension | Definition | Units |
|----------|-----------|------------|-------|
| `FREQ(NFREQ+1)` | 101 | Bin boundaries | h*nu/(k*UTP) |
| `FREQMN(NFREQ)` | 100 | Bin centres (mean frequency) | h*nu/(k*UTP) |
| `WEIGHT(130)` | 100 | bandwidth * FREQMN^3 | (h*nu/kT)^4 |

**Physical frequency:** `nu = FREQMN(L) * UFREQ` where `UFREQ = BOLTZK*UTP/(2*PI*HPLANC) ≈ 2.084e15 Hz`

### Radiation Moments in Y Array (Time-Evolved)
The solution vector `Y` stores radiation moments as **occupation numbers** (dimensionless):

| Index Formula | Quantity | Frame |
|---------------|----------|-------|
| `Y(NVARS*NZON + (NZON-NCND)*(L-1) - NCND + Km)` | J occupation number (zeroth moment) | lab-frame |
| `Y(NVARS*NZON + (NZON-NCND)*(L-1) - NCND + KRAD + Km)` | H occupation number (first moment) | lab-frame |

Where:
- `NVARS = 4`, `KRAD = (NZON-NCND)*NFRUS`
- `NCND` = number of opaque (condensed) inner zones
- Occupation number: `n = c^2/(2h*nu^3) * I_nu`

**Physical J_nu:** `J_nu [erg/cm^2/s/sr/Hz] = (2h*UFREQ^3/c^2) * FREQMN(L)^3 * Y_J(K,L)`

### Radiation Moments from Formal Solution (Feau/eddi.trf)
The `Feau` subroutine computes more accurate radiation moments via ray-by-ray formal solution with Doppler shifts:

| Variable | Dimension | Definition | Frame |
|----------|-----------|------------|-------|
| `FradJ(K,Lfr)` | (Mzon, Nfreq) | J_nu (zeroth moment) | comoving |
| `FradH(K,Lfr)` | (Mzon, Nfreq) | H_nu (first moment) | comoving |
| `FradK(K)` | (Mzon) | K_nu (second moment, per-frequency) | comoving |
| `EddN(K,Lfr)` | (Mzon, Nfreq) | Eddington factor f = K/J | dimensionless |
| `HEdN(Nfreq)` | (Nfreq) | H/J at outer boundary | dimensionless |
| `EDDJ(Mzon,Nfreq)` | (Mzon, Nfreq) | Older Eddington factor array | dimensionless |

**Key:** `EddN(K,Lfr) = FradK(K) / FJLfr(K)` (eddi.trf line 657) — exactly f = K/J.

### Frame Considerations
- **Y array moments**: lab-frame (Eulerian), evolved by moment equations in traneq.f
- **FradJ/FradH/EddN**: comoving-frame, from formal solution with Doppler-shifted frequency interpolation
- **CMFGEN EDDFACTOR**: expects comoving-frame J_nu

**Recommendation:** Use FradJ (comoving) for EDDFACTOR. If FradJ is stale (Feau not called recently), fall back to Y-array J with lab→comoving transformation.

### Bolometric Quantities Already Computed
From lbalsw.f lines 894-931:
```
SUMJ = sum(Y_J(K,L) * WEIGHT(L))  →  TpRAD = (15/pi^4 * SUMJ)^0.25
Flum = CLUMF * sum(Y_H(K,L) * WEIGHT(L)) * R^2  →  luminosity
ERAD = RadC * SUMJ / rho  →  radiation energy density per unit mass
```

---

## 3. CMFGEN EDDFACTOR File Format

### Source Files
- `/ptmp/sbrennan/SN2024pba/cmfgen_source/new_main/subs/open_rw_eddfactor.f`
- `/ptmp/sbrennan/SN2024pba/cmfgen_source/new_main/mod_subs/comp_j_blank.f`
- `/ptmp/sbrennan/SN2024pba/cmfgen_source/new_main/mod_subs/eddfac_rec_defs_mod.f`

### File Format: Direct-Access Binary

| Record | Name | Content |
|--------|------|---------|
| 1 | `DIE_CONT_REC` | Access pointer (integer) |
| 2 | `RV_REC` | R,V coordinate pointer |
| 3 | `EDD_CONT_REC` | Continuum access pointer + NCF + NDEXT |
| 4 | `LINE_CONT_REC` | Line access pointer |
| 5 | `FINISH_REC` | Nonzero when all frequencies computed |
| 6-9 | reserved | — |
| 10-14 | `INITIAL_RV_REC` | R(1:ND_EXT), V(1:ND_EXT), LANG_COORD(1:ND_EXT) |
| 15+ | `INITIAL_ACCESS_REC` | Per-frequency: RJ(1:ND_EXT) + FL (frequency in Hz) |

**Record size:** `WORD_SIZE * (ND_EXT + 1) / UNIT_SIZE` (in direct-access units)

**Companion file:** `EDDFACTOR_INFO` contains ND_EXT, record_size, date.

### What EDDFACTOR Stores
Despite the name, EDDFACTOR stores **J_nu (mean intensity)**, not Eddington factors directly.
Each per-frequency record contains:
- `RJ(I)` for I=1,ND_EXT (J at each depth point, **outer-to-inner**)
- `FL` (frequency in Hz, ordered **highest to lowest**)

### How CMFGEN Uses EDDFACTOR

#### When `COMPUTE_EDDFAC = .TRUE.` (VADAT: `COMP_F=T`, default)
- RJ initialised to 0
- J and Eddington factor f computed self-consistently via `JFEAU_IBC`
- Results written to EDDFACTOR for reuse in subsequent iterations

#### When `COMPUTE_EDDFAC = .FALSE.` (VADAT: `COMP_F=F`)
- RJ read from EDDFACTOR file
- Used as **initial guess** for J iteration
- Eddington factor f computed from the read RJ
- **This is the mode for providing STELLA-derived J**

#### When `USE_FIXED_J = .TRUE.` (VADAT: `USE_FIXED_J=T`)
- RJ read from EDDFACTOR and **kept fixed** (no iteration)
- `K_MOM = RJ/3` (isotropic assumption)
- Requires `RD_LAMBDA = .TRUE.`
- **Not recommended for production** — too restrictive

### Current VADAT Configuration (from VADAT_INIT)
```
T  [COMP_F]       ! Compute Eddington factors at start (default)
F  [USE_FIXED_J]  ! Don't use fixed J
T  [JC_W_EDD]     ! Compute continuum with Eddington factors (opens EDDFACTOR file)
```

### CMFGEN J_nu Units
RJ is in the same units as the source function `TC = eta/chi`:
- For LTE: `TC = B_nu(T) = (2h*nu^3/c^2) / (exp(h*nu/kT) - 1)` [erg/cm^2/s/sr/Hz]
- So **RJ is in erg/cm^2/s/sr/Hz**

---

## 4. STELLA → CMFGEN Variable Mapping

### Physical Conversion

| STELLA Quantity | Conversion | CMFGEN Quantity | Units |
|----------------|------------|-----------------|-------|
| `Y_J(K,L)` (occupation number) | `J_nu = (2h*UFREQ^3/c^2) * FREQMN(L)^3 * Y_J` | `RJ(I)` | erg/cm^2/s/sr/Hz |
| `FREQMN(L)` (dimensionless freq) | `nu = FREQMN(L) * UFREQ` | `FL` | Hz |
| `FradJ(K,L)` (comoving J occ. number) | `J_nu = (2h*UFREQ^3/c^2) * FREQMN(L)^3 * FradJ` | `RJ(I)` | erg/cm^2/s/sr/Hz |
| `EddN(K,L)` (dimensionless) | Direct (no conversion needed) | Eddington factor f | dimensionless |

### Key Constants
```
UFREQ = BOLTZK * UTP / (2*PI*HPLANC) = k * 1e5 / h ≈ 2.084e15 Hz
J_scale = 2h * UFREQ^3 / c^2 ≈ 0.1335  [erg/cm^2/s/sr/Hz per occ.number * FREQMN^3]
```

### Grid Mapping
- **STELLA**: zones 1→NZON (inner→outer), NFREQ frequency groups
- **CMFGEN**: depths 1→ND (outer→inner), NCF frequencies (highest→lowest)
- **Reversal needed**: STELLA inner→outer must be reversed to CMFGEN outer→inner
- **Frequency interpolation**: STELLA's 100 groups → CMFGEN's frequency grid
- **Depth interpolation**: STELLA's zones → CMFGEN's depth grid (use radius or optical depth)

### Frame Considerations
- STELLA `FradJ` is **comoving-frame** (from formal solution with Doppler shifts)
- STELLA `Y_J` is **lab-frame** (from moment equations)
- CMFGEN expects **comoving-frame** J_nu
- **Prefer FradJ** for EDDFACTOR construction

---

## 5. Implementation Plan

### Step 1: Modify STELLA .swd Writer
Add a separate file (e.g., `.swd_rad`) written alongside `.swd` containing:
- Version header
- Frequency grid: FREQ, FREQMN, WEIGHT, UFREQ
- Per-zone per-frequency: FradJ, FradH, EddN (from last Feau call)
- Per-zone per-frequency: Y_J, Y_H (always current, lab-frame)
- Conversion constants: UTP, UR, UTIME, etc.

**Rationale for separate file:** Avoid breaking existing .swd readers.

### Step 2: Build Converter (`convert_stella_swd_to_edd_factor.py`)
1. Read `.swd` for spatial structure (R, V, T, rho)
2. Read `.swd_rad` for frequency-dependent radiation
3. Convert occupation numbers → physical J_nu
4. Map STELLA depth grid → CMFGEN depth grid (radius-based interpolation)
5. Interpolate/rebin STELLA frequency groups → CMFGEN frequency grid
6. Reverse ordering (inner→outer to outer→inner)
7. Write EDDFACTOR binary file + EDDFACTOR_INFO

### Step 3: VADAT Configuration
Set `COMP_F=F` to read pre-computed J from EDDFACTOR instead of computing from scratch.

### Step 4: Test with 1-Day STELLA Model
- Create minimal config
- Run STELLA ~1 day
- Generate EDDFACTOR
- Run CMFGEN control (COMP_F=T) vs experiment (COMP_F=F with STELLA EDDFACTOR)
- Compare convergence

---

## 6. Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| FradJ may be stale (Feau not called at .swd output time) | Also write Y_J; converter can choose |
| Frequency grid mismatch | Explicit interpolation with documentation |
| Depth grid mismatch | Radius-based interpolation, not index-based |
| Non-monotonic velocity | Use radius (always monotonic) for depth mapping |
| Unit conversion errors | Validate against bolometric quantities (T_rad, Flum) |
| EDDFACTOR binary format incompatibility | Test with actual CMFGEN run |
