#!/usr/bin/env python3
"""
NUMERICAL scout of the open Prop `ArchPosDef` (Flint / Discovery 3, 2026-11-12).
Labels: NUMERICAL only. A finite-window numerical failure KILLS the Prop as
stated (it is a universal statement over finite configurations); a finite-window
success would NOT prove it and would NOT prove RH.

VERBATIM SOURCE (GPPVerify @ 5fb6c8d0, branch heart/archpd, read-only clone
~/workspace/gppverify-archpd-readonly):

GppVerify/RiemannHypothesis/WeilGramArchPD.lean:36-37
  noncomputable def archQ {n : Nat} (x c : Fin n -> R) : R := kform archTerm x c
GppVerify/RiemannHypothesis/WeilGramArchPD.lean:44-46
  def ArchPosDef (L : R) : Prop :=
    forall (n : Nat) (x c : Fin n -> R), (forall i j, |x i - x j| <= L) ->
      sum i, c i * Real.cosh (x i / 2) = 0 -> c != 0 -> 0 < archQ x c
GppVerify/RiemannHypothesis/WeilGramArchPD.lean:49
  def ArchPosDefAll : Prop := forall L > 0, ArchPosDef L

GppVerify/RiemannHypothesis/KreinBridge.lean:53-54
  def kform (K : R -> R) {n : Nat} (x : Fin n -> R) (c : Fin n -> R) : R :=
    sum i, sum j, c i * c j * K (x i - x j)

GppVerify/RiemannHypothesis/KreinHPsi.lean:57-59 (archTerm; bigPsi at :46-47,
bigPsi t = if t = 0 then 1 else (sin(t/2)/(t/2))^2)
  noncomputable def archTerm (d : R) : R :=
    (1 / pi) * integral t in Ioi 0, bigPsi t^2 * cos(t*d) *
      ((logDeriv Gamma (1/4 + (t/2)*I)).re - log pi)

So: A_ij = archTerm(x_i - x_j), ell_c(c) = sum_i c_i cosh(x_i/2), and the
window condition is PAIRWISE SPREAD |x_i - x_j| <= L (not x_i in [-L, L]).
archTerm is even, so A_ij depends on |x_i - x_j|.

VALIDATION ANCHOR: HPsiPlan.lean header (same commit) records, for n=2,
x=(0,a), c=(cosh(a/2),-1) in ker ell_c, the mpmath values of
kform(poleTerm+archTerm) = -0.0054, -0.205, -2.00, -6.85, -49.7 at
a = 0.05, 0.3, 1, 2, 4.  This script recomputes pole+arch from its own
archTerm and checks agreement; the pole part on ker ell_c is the PROVED
Lean identity kform_poleTerm_of_ker = -2*Psi(i/2)^2*(sum c_i sinh(x_i/2))^2.

Method: scan with scipy (float, complex digamma) at dps-equivalent ~1e-11;
confirm the n=2 falsifier series and one representative per L with mpmath
at 35-50 digits, including the exact Rayleigh quotient of the falsifying
vector.  Eigenvalue error from entry error eps is <= n*eps (Weyl).
"""
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import digamma

LOGPI = np.log(np.pi)

def psi_big(t):
    if t == 0.0:
        return 1.0
    return (np.sin(t / 2) / (t / 2)) ** 2

def arch_scipy(d, T=250.0):
    d = abs(float(d))
    def f(t):
        if t == 0.0:
            w = digamma(0.25).real - LOGPI
            return w  # Psi(0)^2 * cos(0) * w
        return psi_big(t) ** 2 * np.cos(t * d) * (digamma(0.25 + 0.5j * t).real - LOGPI)
    val, err = quad(f, 0, T, limit=400, epsabs=1e-12, epsrel=1e-12,
                    points=[1, 2, 5, 10, 20, 40, 80, 160])
    # crude tail bound: Psi(t)^2 <= 16/t^4, |Re psi(1/4+it/2)-log pi| <= log t (t>=250, checked vs samples)
    tail = (16 / np.pi) * (np.log(T) / (3 * T ** 3) + 1 / (9 * T ** 3))
    return val / np.pi, err / np.pi + tail

_cache = {}
def A_mat(x):
    x = np.asarray(x, float)
    n = len(x)
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            key = round(abs(x[i] - x[j]), 10)
            if key not in _cache:
                _cache[key] = arch_scipy(key)[0]
            A[i, j] = _cache[key]
    return A

def proj_eigmin(x):
    """smallest eigenvalue of A restricted to ker ell_c, and full (off-ker) eigmin."""
    x = np.asarray(x, float)
    n = len(x)
    A = A_mat(x)
    ell = np.cosh(x / 2)
    # orthonormal basis of ker ell via SVD
    _, _, Vt = np.linalg.svd(ell.reshape(1, -1))
    Z = Vt[1:].T  # n x (n-1)
    M = Z.T @ A @ Z
    w = np.linalg.eigvalsh(M)
    wf = np.linalg.eigvalsh(A)
    # minimizing vector in ker
    evals, evecs = np.linalg.eigh(M)
    c = Z @ evecs[:, 0]
    return w[0], wf[0], c

def configs(L, n):
    half = (n - 1) / 2
    equi = np.linspace(-L / 2, L / 2, n)
    clust = (np.arange(n) - half) * (0.1 * L / max(n - 1, 1))
    ends = np.array([-L / 2 if i % 2 == 0 else L / 2 for i in range(n)])
    return {"equispaced": equi, "clustered0.1L": clust, "endpoints": ends}

def main():
    out = {"label": "NUMERICAL", "results": [], "n2": []}
    print("=== ArchPosDef scout (NUMERICAL) ===")
    print("archTerm(0) scipy =", arch_scipy(0.0))
    # ---- n=2 smallest falsifier series + validation anchor ----
    bigPsiI2 = 16 * np.sinh(0.25) ** 2
    anchor = {0.05: -0.0054, 0.3: -0.205, 1.0: -2.00, 2.0: -6.85, 4.0: -49.7}
    print("\n-- n=2, x=(0,a), c=(cosh(a/2),-1) in ker ell_c --")
    for a, ref in anchor.items():
        A = A_mat([0.0, a])
        c = np.array([np.cosh(a / 2), -1.0])
        archQ = float(c @ A @ c)
        pole = -2 * bigPsiI2 ** 2 * (c[1] * np.sinh(a / 2)) ** 2
        ray = archQ / float(c @ c)
        print(f"a={a:>5}: archQ={archQ:+.10f}  Rayleigh={ray:+.10f}  "
              f"pole+arch={pole + archQ:+.6f} (HPsiPlan anchor {ref})")
        out["n2"].append({"a": a, "archQ": archQ, "rayleigh": ray,
                          "pole_plus_arch": pole + archQ, "anchor": ref})
    # ---- scan n=3..8 ----
    print("\n-- scan: projected (ker ell_c) smallest eigenvalue / full-matrix smallest eigenvalue --")
    hdr = f"{'L':>4} {'config':>14} {'n':>2} {'diam':>6} {'lam_min_ker':>14} {'lam_min_full':>14} {'dimker':>6}"
    print(hdr)
    for L in [0.5, 1.0, 2.0, 4.0]:
        for n in range(3, 9):
            for name, x in configs(L, n).items():
                lam, lamf, c = proj_eigmin(x)
                diam = float(np.max(x) - np.min(x))
                print(f"{L:>4} {name:>14} {n:>2} {diam:>6.2f} {lam:>14.8f} {lamf:>14.8f} {n - 1:>6}")
                out["results"].append({"L": L, "config": name, "n": n, "diam": diam,
                                        "lam_min_ker": float(lam), "lam_min_full": float(lamf),
                                        "dim_ker": n - 1})
    # ---- mpmath high-precision confirmation ----
    print("\n-- mpmath confirmation (dps=40) --")
    import mpmath as mp
    mp.mp.dps = 40
    def arch_mp(d):
        d = mp.mpf(d)
        f = lambda t: (mp.mpf(1) if t == 0 else (mp.sin(t / 2) / (t / 2)) ** 2) ** 2 \
            * mp.cos(t * d) * (mp.re(mp.digamma(mp.mpc(mp.mpf("0.25"), t / 2))) - mp.log(mp.pi))
        return (1 / mp.pi) * mp.quad(f, [0, mp.inf])
    conf = []
    A0 = arch_mp(0)
    print("archTerm(0) mp dps40 =", A0)
    for a in ["0.3", "1", "2"]:
        Aa = arch_mp(a)
        av = mp.mpf(a)
        c0, c1 = mp.cosh(av / 2), mp.mpf(-1)
        q = c0 * c0 * A0 + 2 * c0 * c1 * Aa + c1 * c1 * A0
        ray = q / (c0 * c0 + 1)
        print(f"a={a}: archQ={mp.nstr(q, 20)}  Rayleigh={mp.nstr(ray, 20)}")
        conf.append({"a": a, "archQ_mp": mp.nstr(q, 20), "rayleigh_mp": mp.nstr(ray, 20)})
    # representative projected case: L=1, equispaced n=4, mpmath Rayleigh of scipy minimizer
    L, n = 1.0, 4
    x = np.linspace(-L / 2, L / 2, n)
    lam, lamf, cvec = proj_eigmin(x)
    lags = sorted({round(abs(float(x[i] - x[j])), 10) for i in range(n) for j in range(n)})
    Amp = {lag: arch_mp(lag) for lag in lags}
    cm = [mp.mpf(float(v)) for v in cvec]
    # exact projection onto ker ell_c in mpmath arithmetic
    ellm = [mp.cosh(mp.mpf(float(v)) / 2) for v in x]
    dot = sum(cm[i] * ellm[i] for i in range(n)); nn = sum(e * e for e in ellm)
    cm = [cm[i] - ellm[i] * dot / nn for i in range(n)]
    q = sum(cm[i] * cm[j] * Amp[round(abs(float(x[i] - x[j])), 10)] for i in range(n) for j in range(n))
    nrm = sum(v * v for v in cm)
    print(f"L=1 equispaced n=4: scipy lam_min_ker={lam:+.10f}; mpmath Rayleigh of (exactly projected) "
          f"minimizer = {mp.nstr(q / nrm, 20)}; ker residual={mp.nstr(sum(cm[i]*ellm[i] for i in range(n)), 5)}")
    conf.append({"case": "L=1 equispaced n=4", "scipy_lam": lam, "mp_rayleigh": mp.nstr(q / nrm, 20)})
    out["mpmath_confirmation"] = conf
    print("\nJSON:", json.dumps(out)[:200], "...")
    return out

if __name__ == "__main__":
    main()
