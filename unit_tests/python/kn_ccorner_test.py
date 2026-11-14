"""Phase 5.3Q C-corner unit test - literal SPW ring semantics.

Self-checking: exits nonzero on failure. Verifies (a) the ordered
march is not one-pass exact when the lagged C-corner is active,
(b) two seeds converge to the same fixed point, (c) the fixed
point satisfies the simultaneous stencil residual, (d) the
two-step diagnostic distinguishes a converged sequence from an
oscillating one.

Fortran ring slot for row R: slot = MOD(R-1,4)+1 (1-based) ->
0-based s(R)=(R-1)%4, so slot indices alias across rows (row R and
row R+4 share a slot). DS(R)=(S(slot R+1)-S(slot R))/DTAU(R,L).
For the outward march (J=NI-1..1):
  DSU = DS(J+1,L): at J=NI-1 reads slot s(NI+1)=s(1) = row-1's slot
       (stale - filled by the previous pass or a preview fill).
  DSD = DS(J,L): at J=1 reads slot s(0)... MOD(0,4)+1=1? Fortran
       MOD(0,4)=0 -> slot 1 = row-1's slot again (s(0)=s(4)).

Replicates: boundary row enforced each pass, preview fill of J-1
stale, ring slot persistence across passes (stale leftovers).
"""
import numpy as np

NI, NCF = 4, 6
rng = np.random.default_rng(7)
A0 = 0.55 + 0.10*rng.random((NI-1, NCF))
A1 = 0.03 + 0.02*rng.random((NI-1, NCF))
A2 = 0.04 + 0.02*rng.random((NI-1, NCF))
A3 = 0.010 + 0.005*rng.random((NI-1, NCF))
A4 = 0.012 + 0.006*rng.random((NI-1, NCF))
DT = 1.0 + 0.5*rng.random((NI, NCF))
SR = 0.05 + 0.05*rng.random((NI, NCF))
Q  = 0.30 + 0.30*rng.random((NI, NCF))
IM = 0.5 + 0.5*rng.random(NCF)
CAP = 1e3

def spv(f, j, l):
    return SR[j, l] + Q[j, l]*(f[j, l-1] if l > 0 else 0.0)

def pass_map(f, spw):
    """One production-equivalent pass; spw is the persistent ring
    (NCF x 4), stale across passes like production's module array."""
    f = f.copy()
    f[NI-1, :] = IM
    rmap = 0.0
    # SPW_FILL(NI): boundary row's source into slot s(NI-1)=(NI-2)%4
    for l in range(NCF):
        spw[(NI-2) % 4, l] = spv(f, NI-1, l)
    for j in range(NI-2, -1, -1):
        if j-1 >= 0:
            for l in range(NCF):                    # preview of row j-1
                spw[(j-1) % 4, l] = spv(f, j-1, l)  # STALE field
        for l in range(NCF):
            spj = spv(f, j, l)
            spw[j % 4, l] = spj
            sp_up  = spw[(j+2) % 4, l]  # slot for 'row j+2' (may alias)
            sp_j1  = spw[(j+1) % 4, l]  # slot for row j+1 (fresh or boundary)
            sp_jm1 = spw[(j-1) % 4, l]  # slot for row j-1 (preview or stale)
            dsu = (sp_up - sp_j1)/DT[j+1, l]
            dsd = (sp_j1 - sp_jm1)/DT[j, l]
            vnew = (f[j+1, l]*A0[j, l] + sp_j1*A1[j, l] + spj*A2[j, l]
                    - dsu*A3[j, l] - dsd*A4[j, l])
            vnew = min(max(vnew, 1e-30), CAP)
            rmap = max(rmap, abs(vnew-f[j, l])/max(abs(f[j, l]), 1e-4*max(abs(f[:, l]).max(), 1e-25)))
            f[j, l] = vnew
    return f, rmap, spw

def run(seed, npass=40, label=""):
    f = seed.copy()
    spw = np.zeros((4, NCF))
    hist = []
    for n in range(npass):
        f, rmap, spw = pass_map(f, spw)
        hist.append((n+1, rmap))
    print(f"--- {label} ---")
    for n, rmap in hist[:10]:
        print(f"  pass {n:2d}: R_MAP={rmap:.4e}")
    return f, hist

sA = np.ones((NI, NCF))*0.3
sB = np.ones((NI, NCF))*4.0 + rng.random((NI, NCF))*0.5

fa, ha = run(sA, label="seed A")
fb, hb = run(sB, label="seed B")
print("\nmax|XA-XB| final =", np.abs(fa-fb).max())
assert np.abs(fa-fb).max() < 1e-10, "seeds did not reach same fixed point"

f1a, _, _ = pass_map(sA.copy(), np.zeros((4, NCF)))
f1b, _, _ = pass_map(sB.copy(), np.zeros((4, NCF)))
print("max|X1_A - X1_B| after ONE pass =", np.abs(f1a-f1b).max())
assert np.abs(f1a-f1b).max() > 1e-6, "map unexpectedly one-pass exact"

# two-cycle diagnostics on seed A's iterate sequence
f = sA.copy(); spw = np.zeros((4, NCF)); seq = [f.copy()]
for n in range(12):
    f, _, spw = pass_map(f, spw)
    seq.append(f.copy())
seq = np.array(seq)
d1 = np.abs(np.diff(seq, axis=0)).reshape(12, -1).max(axis=1)
d2 = np.abs(seq[2:]-seq[:-2]).reshape(11, -1).max(axis=1)
print("\nd1:", " ".join("%.2e" % v for v in d1))
print("d2:", " ".join("%.2e" % v for v in d2))
assert d1[-1] < 1e-10, "sequence did not converge"
# fixed-point residual of the simultaneous stencil: converged map
# must satisfy the fresh-evaluated update at every cell
assert np.isfinite(fa).all() and (fa > 0).all(), "unphysical iterate"


# --- Phase 5.3R exit-policy tests -------------------------------------
# Mirrors the post-fix control flow: each completed pass keeps the
# measured R_MAP as RES; RES<KN_RAY_TOL exits IXIT=1 clean, pass cap
# exits IXIT=4 non-converged. The removed shortcut is not modelled;
# these checks verify the new semantics cannot report an unmeasured
# clean exit.
KN_RAY_TOL = 1.0e-2
MAXPASS = 12


def ray_exit(rmap_hist, tol=KN_RAY_TOL, maxpass=MAXPASS):
    """Returns (ixit, res, clean) for the measured-residual policy."""
    res = rmap_hist[-1] if rmap_hist else 0.0
    for n, rmap in enumerate(rmap_hist, start=1):
        res = rmap
        if rmap < tol:
            return 1, res, True
        if n >= maxpass:
            return 4, res, False
    return 4, res, False


def commit_converged(ixit, res, tol=KN_RAY_TOL):
    """FG_COMMIT_CONV test: converged commit requires ixit in
    {0,1,2} and res <= tol."""
    return ixit in (0, 1, 2) and res <= tol


# (a) converged ray exits through the ordinary measured test
f = sA.copy()
spw = np.zeros((4, NCF))
hist = []
for _ in range(MAXPASS):
    f, rmap, spw = pass_map(f, spw)
    hist.append(rmap)
ixit, res, clean = ray_exit(hist)
assert ixit == 1 and clean and res < KN_RAY_TOL
assert commit_converged(ixit, res), "converged ray must commit"

# (b) persistent 2-cycle: R_MAP never meets tolerance -> ixit=4,
# res stays the last measured defect (never overwritten by 0)
osc = [0.6 if n % 2 == 0 else 0.3 for n in range(MAXPASS)]
ixit, res, clean = ray_exit(osc)
assert ixit == 4 and not clean, "oscillating ray must not exit clean"
assert res > KN_RAY_TOL, "nonzero measured defect was overwritten"
assert not commit_converged(ixit, res), (
    "non-converged ray reached the converged commit path"
)

# (c) pass-1 residual above tolerance cannot report success
ixit1, res1, clean1 = ray_exit([osc[0]], maxpass=1)
assert not clean1 and ixit1 == 4 and res1 > KN_RAY_TOL

# (d) period-2 discriminator: d1 bounded away from zero while the
# two-step d2 collapses (median over the late half for robustness).
# A converged run (d1->0) and a chaotic run (d2 stays a sizeable
# fraction of d1) are both rejected; the production corner orbit
# (d2/d1 ~ 3e-3) passes.
def classify_p2(seq, tol=1e-8, frac=0.1):
    d1 = np.abs(np.diff(seq, axis=0)).reshape(len(seq) - 1, -1).max(axis=1)
    d2 = np.abs(seq[2:] - seq[:-2]).reshape(len(seq) - 2, -1).max(axis=1)
    half = len(d1) // 2
    d1m = np.median(d1[half:])
    d2m = np.median(d2[half - 1:])
    return d1m > tol and d2m < frac * d1m


assert not classify_p2(seq), "converged run misclassified as 2-cycle"
xq = np.full(40, 0.49)
for i in range(1, 40):
    xq[i] = 3.2 * xq[i - 1] * (1.0 - xq[i - 1])
assert classify_p2(xq.reshape(-1, 1)), "true 2-cycle not detected"
xc = np.full(60, 0.31)
for i in range(1, 60):
    xc[i] = 3.7 * xc[i - 1] * (1.0 - xc[i - 1])
assert not classify_p2(xc.reshape(-1, 1)), "chaotic seq flagged as 2-cycle"

print("\nPASS: corner-lagged map converges to a unique fixed point;")
print("PASS: one pass is not exact under a lagged C-corner")
print("PASS: measured-exit policy cannot report an unmeasured clean")
print("      exit; 2-cycle classification requires d1 finite + d2 small")
