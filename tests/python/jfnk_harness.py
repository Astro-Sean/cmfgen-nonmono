#!/usr/bin/env python3
"""Data-free JFNK verification harness.

Compiles the production Fortran JFNK solver (newton_hook_v1 +
gmres_m_v1 in new_main/mod_subs/jfnk_solver_v1.f90) together with a
thin CLI driver via subprocess/ifort, then runs three analytical
tests and produces matplotlib plots:

  T1  finite-difference "U-curve": FD J*v vs exact J*v over eps
  T2  outer Newton convergence on the 2D Rosenbrock system
  T3  inner GMRES convergence on an ill-conditioned system, with and
      without a Jacobi preconditioner

To plug a different Fortran solver in, edit the CONFIG block below
and match the CLI protocol implemented in jfnk_cli_driver.f90
(commands: fdjvp / newton / gmres).  Any interface with the same
stdout contract drops in unchanged.
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# CONFIG: sources and compile line.
# Replace FORTRAN_SOURCES with your own module(s) when plugging in a
# different JFNK implementation; the driver only needs the module
# symbols it USEs (jfnk_solver_v1_mod) plus the stub int_mod.
# ----------------------------------------------------------------------
TREE = Path(__file__).resolve().parents[2]
TESTDIR = Path(__file__).resolve().parents[1]
FORTRAN_SOURCES = [
    TESTDIR / "f90" / "jfnk_int_mod_stub.f90",
    TREE / "new_main" / "mod_subs" / "jfnk_solver_v1.f90",
    TESTDIR / "f90" / "jfnk_cli_driver.f90",
]
BINARY = TESTDIR / "f90" / "bin" / "jfnk_cli"
OUTDIR = TESTDIR / "plots" / "harness"

IFORT = shutil.which("ifort")
GFORT = shutil.which("gfortran")
if IFORT:
    COMPILE = [
        IFORT, "-qopenmp", "-mkl", "-O0", "-g", "-traceback",
        "-diag-disable=10448,10412,10434",
        "-module", str(BINARY.parent),
    ]
elif GFORT:
    COMPILE = [GFORT, "-fopenmp", "-O0", "-g", "-fbacktrace",
               "-llapack", "-lblas", f"-J{BINARY.parent}"]
else:
    sys.exit("no Fortran compiler found (need ifort or gfortran)")


def ensure_binary():
    """Compile the solver + CLI driver once (cached by mtime)."""
    if BINARY.exists():
        newest = max(p.stat().st_mtime for p in FORTRAN_SOURCES)
        if BINARY.stat().st_mtime > newest:
            return
    BINARY.parent.mkdir(parents=True, exist_ok=True)
    cmd = COMPILE + [str(p) for p in FORTRAN_SOURCES] + [
        "-o", str(BINARY)]
    print("compile:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_cli(*args):
    """Run the CLI driver, return stdout."""
    r = subprocess.run(
        [str(BINARY)] + [str(a) for a in args],
        capture_output=True, text=True, check=True)
    return r.stdout


# ----------------------------------------------------------------------
# Test 1: finite-difference U-curve
# ----------------------------------------------------------------------
def rosen_f(x):
    return np.array([1.0 - x[0], 100.0 * (x[1] - x[0] ** 2)])


def rosen_jac(x):
    return np.array([[-1.0, 0.0], [-200.0 * x[0], 100.0]])


def test_fd_ucurve():
    x = np.array([-1.2, 1.0])
    v = np.array([1.0, 1.0]) / np.sqrt(2.0)
    jv_exact = rosen_jac(x) @ v
    epss = np.logspace(-2, -14, 49)
    errs = []
    for eps in epss:
        out = run_cli("fdjvp", f"{eps:.16e}", *x, *v)
        jv_fd = np.array([float(t) for t in out.split()[1:3]])
        errs.append(
            np.linalg.norm(jv_fd - jv_exact) / np.linalg.norm(jv_exact))
    errs = np.array(errs)
    imin = np.argmin(errs)
    passed = errs[imin] < 1e-6 and 1e-9 <= epss[imin] <= 1e-7

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(epss, errs, "o-", ms=4)
    ax.axvline(np.sqrt(np.finfo(float).eps), color="r", ls="--",
               label=r"$\sqrt{\epsilon_{mach}}$")
    ax.set_xlabel(r"finite-difference step scale $\epsilon$")
    ax.set_ylabel(r"rel. error in FD $J v$")
    ax.set_title("T1: FD Jacobian-vector product U-curve "
                 "(production step rule)")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    out = OUTDIR / "t1_fd_ucurve.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  min err {errs[imin]:.2e} at eps {epss[imin]:.2e} -> {out}")
    return passed, "min err < 1e-6 and eps* in [1e-9, 1e-7]"


# ----------------------------------------------------------------------
# Test 2: outer Newton convergence on Rosenbrock
# ----------------------------------------------------------------------
def test_newton_rosen():
    out = run_cli("newton", -1.2, 1.0)
    info = int(out.split("INFO")[1].split()[0])
    xline = [ln for ln in out.splitlines() if ln.startswith("X ")][0]
    x = np.array([float(t) for t in xline.split()[1:]])
    hist = np.array([
        [float(t) for t in ln.split()]
        for ln in out.split("\nHIST\n")[1].strip().splitlines()])
    its, res = hist[:, 0], hist[:, 1]
    passed = (info == 0 and np.allclose(x, [1.0, 1.0], atol=1e-8)
              and res[-1] < 1e-10 and len(res) <= 15)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(its, res, "o-", ms=5)
    ax.axhline(1e-15, color="r", ls="--", label="machine floor")
    ax.set_xlabel("outer Newton iterate")
    ax.set_ylabel(r"$\|F(x)\|_2$")
    ax.set_title("T2: newton_hook_v1 on Rosenbrock "
                 f"(x0=(-1.2,1), info={info})")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    out = OUTDIR / "t2_newton_rosen.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  converged x={x}, |F|={res[-1]:.2e} in "
          f"{int(its[-1])} iterates -> {out}")
    return passed, "info=0, x=(1,1), |F|<1e-10, <=15 iterates"


# ----------------------------------------------------------------------
# Test 3: inner GMRES with/without Jacobi preconditioner
# ----------------------------------------------------------------------
def test_gmres_pc():
    rng = np.random.default_rng(20261008)
    n, m, itmax = 100, 20, 150
    # Ill-scaled system: A = D*B, D log-uniform over 8 decades and
    # B = I + 0.3*S well-conditioned (cond ~2).  Jacobi recovers B
    # exactly under left preconditioning, so the PC Krylov problem is
    # easy while unpreconditioned A has cond ~1e8.
    s = rng.uniform(-1.0, 1.0, (n, n))
    s = 0.5 * (s + s.T)
    s = s / np.linalg.norm(s, 2)
    b_mat = np.eye(n) + 0.3 * s
    d = 10.0 ** rng.uniform(0.0, 8.0, n)
    a = np.diag(d) @ b_mat
    x_true = rng.uniform(-1.0, 1.0, n)
    b = a @ x_true
    print(f"  n={n} cond(A)~{np.linalg.cond(a):.2e} "
          f"cond(M^-1 A)~{np.linalg.cond(b_mat):.1f}")

    curves = {}
    for label, pc in (("no preconditioner", 0),
                      ("Jacobi (diag) PC", 1)):
        with tempfile.TemporaryDirectory() as td:
            # Fortran list-directed read fills a(:,:) column-major, so
            # write A^T row-major for the driver to see A itself.
            np.savetxt(Path(td) / "a.txt", a.T)
            np.savetxt(Path(td) / "b.txt", b)
            out = run_cli("gmres", n, m, itmax, pc,
                          Path(td) / "a.txt", Path(td) / "b.txt")
        rows = [ln.split() for ln in out.splitlines()
                if ln and ln[0].isdigit()]
        ks = np.array([int(r[1]) for r in rows])
        ls_res = np.array([float(r[2]) for r in rows])
        infos = np.array([int(r[4]) for r in rows])
        conv = np.argmax(infos == 0)
        if infos[conv] == 0:
            ks, ls_res = ks[:conv + 1], ls_res[:conv + 1]
        # relres in the norm each method minimizes (GMRES convention):
        # |Ax-b|/|b| unpreconditioned, |M^-1(Ax-b)|/|M^-1 b| with PC.
        curves[label] = ks, ls_res / ls_res[0]

    fig, ax = plt.subplots(figsize=(7, 5))
    for label, (ks, relres) in curves.items():
        ax.semilogy(ks, relres, ".-", ms=5, label=label)
    ax.set_xlabel("Krylov iteration")
    ax.set_ylabel(r"relative residual (each method's minimized norm)")
    ax.set_title("T3: gmres_m_v1 on ill-scaled A=D*B (left Jacobi PC)")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    out = OUTDIR / "t3_gmres_pc.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)

    it0 = curves["no preconditioner"][0][-1]
    it1 = curves["Jacobi (diag) PC"][0][-1]
    r1_end = curves["Jacobi (diag) PC"][1][-1]
    r0_end = curves["no preconditioner"][1][-1]
    passed = it1 < it0 and r1_end < 1e-8
    print(f"  Jacobi: {it1} its to relres {r1_end:.1e}; "
          f"noPC: {it0} its at relres {r0_end:.1e} -> {out}")
    return passed, "Jacobi converges (<1e-8) in far fewer its"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    ensure_binary()
    results = []
    for name, fn in (("T1 FD U-curve", test_fd_ucurve),
                     ("T2 Newton Rosenbrock", test_newton_rosen),
                     ("T3 GMRES Jacobi PC", test_gmres_pc)):
        print(f"[{name}]")
        try:
            ok, crit = fn()
        except subprocess.CalledProcessError as e:
            print(e.stdout, e.stderr)
            ok, crit = False, "driver crashed"
        print(f"  {'PASS' if ok else 'FAIL'}  ({crit})\n")
        results.append((name, ok))
    print("=" * 50)
    for name, ok in results:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    sys.exit(0 if all(ok for _, ok in results) else 1)


if __name__ == "__main__":
    main()
