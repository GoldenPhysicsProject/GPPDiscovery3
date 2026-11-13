#!/usr/bin/env python3
"""
NUMERICAL scout of the open Prop `SemilocalArchDomAll` (Flint / Discovery 3, 2026-11-12).
Labels: NUMERICAL only. A negative D at any tested instance would KILL SemilocalArchDom
(and hence the banked equivalence would kill H_psi) at that instance; a finite-window
numerical success does NOT prove SemilocalArchDomAll, H_psi, or RH.

=== STEP 1: EXACT STATEMENTS, VERBATIM (GPPVerify @ 5fb6c8d0, heart/archpd,
read-only clone ~/workspace/gppverify-archpd-readonly) ===

WeilGramContraction.lean:85-90 (SemilocalArchDom):
  def SemilocalArchDom (S : Finset Nat) (L : R) : Prop :=
    forall (n : Nat) (x c : Fin n -> R), (forall i j, |x i - x j| <= L) ->
      sum i, c i * Real.cosh (x i / 2) = 0 ->
        kform (primeTermS S) x c + 2 * bigPsiI2 ^ 2 * (sum i, c i * Real.sinh (x i / 2)) ^ 2
          <= kform archTerm x c
WeilGramContraction.lean:125 (SemilocalArchDomAll):
  def SemilocalArchDomAll : Prop := forall L > 0, SemilocalArchDom (primesBelowExp L) L
WeilGramContraction.lean: PROVED hPsi_iff_semilocalArchDomAll : HPsi <-> SemilocalArchDomAll.
WeilGramArchPD.lean:40-42 (primePoleQ) and :36-37 (archQ):
  noncomputable def primePoleQ (S : Finset Nat) {n : Nat} (x c : Fin n -> R) : R :=
    kform (primeTermS S) x c + 2 * bigPsiI2 ^ 2 * (sum i, c i * Real.sinh (x i / 2)) ^ 2
  noncomputable def archQ {n : Nat} (x c : Fin n -> R) : R := kform archTerm x c
  PROVED archDom_iff_le: SemilocalArchDom S L <-> forall ..., primePoleQ S x c <= archQ x c.
WeilGramSemilocal.lean:34-36 (primeTermS):
  noncomputable def primeTermS (S : Finset Nat) (d : R) : R :=
    sum' n : Nat, if n.minFac in S then vonMangoldt n / Real.sqrt n *
      (triG (Real.log n + d) + triG (-Real.log n + d)) else 0
WeilGramSemilocal.lean:46-47 (primesBelowExp):
  noncomputable def primesBelowExp (L : R) : Finset Nat :=
    (Finset.range (ceil (Real.exp (2 + L)))).filter Nat.Prime
  PROVED primeTermS_eq: if Covers S L and |d| <= L then primeTermS S d = primeTerm d,
  and covers_primesBelowExp: Covers (primesBelowExp L) L.
KreinHPsi.lean:50-61: triG supported in [-2,2] (piecewise cubic below), primeTerm the
  unrestricted tsum (finite: only n < e^{2+|d|} contribute), Kpsi = pole - prime + arch.

ORIENTATION: P(c) := primePoleQ S x c is the LEFT (dominated) side, A(c) := archQ x c the
RIGHT side: ArchDom is P(c) <= A(c). This script tests D(c) := A(c) - P(c) >= 0 on ker l_c,
with S = primesBelowExp L (matched window). NOTE: ArchDom does NOT assume ArchPosDef and
does NOT exclude c = 0 (equality holds trivially there); A(c) itself is typically NEGATIVE
(ArchPosDef is KILLED, see flint_archposdef_scout_2026_11_12) -- ArchDom is a comparison
of two negative quantities, not positivity of either.
IDENTITY (PROVED in the files, kform_KS_of_ker + kform_poleTerm_of_ker): on ker l_c,
kform Kpsi x c = A(c) - P(c) = D(c) exactly, once S covers L. The n=2 table cross-checks
D against a direct kform(poleTerm - primeTerm + archTerm) evaluation (independent route).

No supremum/limit over S and no abstract S-matrix block occur: primeTermS is the concrete
finite kernel above, so a faithful finite instantiation exists and is what is computed.
primeTermS with the matched S is evaluated EXPLICITLY (minFac filter) and checked against
the full primeTerm at every lag used (they must agree by primeTermS_eq).
"""
import json
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import digamma

LOGPI = np.log(np.pi)
BIGPSI_I2 = 16 * np.sinh(0.25) ** 2

# ---------- archTerm (scipy, as in the ArchPosDef scout; validated vs HPsiPlan anchors) ----------
def psi_big(t):
    if t == 0.0:
        return 1.0
    return (np.sin(t / 2) / (t / 2)) ** 2

def arch_scipy(d, T=250.0):
    d = abs(float(d))
    def f(t):
        if t == 0.0:
            return digamma(0.25).real - LOGPI
        return psi_big(t) ** 2 * np.cos(t * d) * (digamma(0.25 + 0.5j * t).real - LOGPI)
    val, err = quad(f, 0, T, limit=400, epsabs=1e-12, epsrel=1e-12,
                    points=[1, 2, 5, 10, 20, 40, 80, 160])
    tail = (16 / np.pi) * (np.log(T) / (3 * T ** 3) + 1 / (9 * T ** 3))
    return val / np.pi, err / np.pi + tail

# ---------- prime kernels (exact finite sums) ----------
def triG(y):
    ay = abs(y)
    if ay <= 1:
        return 2 / 3 - y * y + ay ** 3 / 2
    if ay <= 2:
        return (2 - ay) ** 3 / 6
    return 0.0

def _factor_prime(n):
    """Return p if n = p^k (k>=1), else None. Also minFac."""
    m, d, facs = n, 2, []
    while d * d <= m:
        while m % d == 0:
            facs.append(d); m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        facs.append(m)
    if not facs:
        return None, n
    p = facs[0]
    return (p if all(q == p for q in facs) else None), p

def prime_kernel(d, S=None):
    """primeTerm(d); if S is a set of primes, primeTermS S d (n kept iff minFac(n) in S)."""
    d = float(d)
    nmax = int(math.exp(2 + abs(d))) + 2
    s = 0.0
    for n in range(2, nmax + 1):
        p, minfac = _factor_prime(n)
        if p is None:
            continue
        if S is not None and minfac not in S:
            continue
        ln = math.log(n)
        s += math.log(p) / math.sqrt(n) * (triG(ln + d) + triG(-ln + d))
    return s

def primes_below_exp_pure(L):
    N = int(math.ceil(math.exp(2 + L)))
    sieve = [True] * N
    out = set()
    for i in range(2, N):
        if sieve[i]:
            out.add(i)
            for j in range(i * i, N, i):
                sieve[j] = False
    return out

# ---------- matrices ----------
_arch_cache, _prime_cache, _primeS_cache = {}, {}, {}
def lag_val(cache, fn, key):
    if key not in cache:
        cache[key] = fn(key)
    return cache[key]

def build_mats(x, L):
    x = np.asarray(x, float); n = len(x)
    S = primes_below_exp_pure(L)
    A = np.zeros((n, n)); Pm = np.zeros((n, n)); PmS = np.zeros((n, n)); Pole = np.zeros((n, n))
    maxdiff = 0.0
    for i in range(n):
        for j in range(n):
            d = round(abs(x[i] - x[j]), 10)
            A[i, j] = lag_val(_arch_cache, lambda k: arch_scipy(k)[0], d)
            Pm[i, j] = lag_val(_prime_cache, lambda k: prime_kernel(k), d)
            PmS[i, j] = lag_val(_primeS_cache, lambda k: prime_kernel(k, S), d)
            maxdiff = max(maxdiff, abs(Pm[i, j] - PmS[i, j]))
            Pole[i, j] = 2 * BIGPSI_I2 ** 2 * np.cosh((x[i] - x[j]) / 2)
    s = np.sinh(x / 2)
    D = A - PmS - 2 * BIGPSI_I2 ** 2 * np.outer(s, s)   # A - P, P = primePoleQ
    K = Pole - Pm + A                                    # direct Kpsi form matrix
    return A, PmS, Pm, D, K, maxdiff

def proj_eigmin(M, x):
    x = np.asarray(x, float)
    ell = np.cosh(x / 2)
    _, _, Vt = np.linalg.svd(ell.reshape(1, -1))
    Z = Vt[1:].T
    Mp = Z.T @ M @ Z
    w, v = np.linalg.eigh(Mp)
    return w[0], w[-1], Z @ v[:, 0]

def configs(L, n):
    half = (n - 1) / 2
    return {"equispaced": np.linspace(-L / 2, L / 2, n),
            "clustered0.1L": (np.arange(n) - half) * (0.1 * L / max(n - 1, 1)),
            "endpoints": np.array([-L / 2 if i % 2 == 0 else L / 2 for i in range(n)])}

def main():
    out = {"label": "NUMERICAL", "n2": [], "scan": [], "primeS_vs_prime_maxdiff": 0.0}
    print("=== SemilocalArchDom scout (NUMERICAL) — D(c) = archQ - primePoleQ on ker l_c, S = primesBelowExp L ===")
    print(f"bigPsiI2 = {BIGPSI_I2:.12f}; archTerm(0) scipy = {arch_scipy(0.0)}")
    # ---- n=2 smallest falsifier family (the family that killed ArchPosDef) ----
    print("\n-- n=2, x=(0,a), c=(cosh(a/2),-1) in ker l_c, L=a (matched) --")
    print(f"{'a':>5} {'archQ':>13} {'primeForm':>13} {'poleCorr':>10} {'P=primePoleQ':>13} {'D=A-P':>11} {'Kpsi_direct':>12} {'D/|c|^2':>10} {'D/|A|':>9}")
    for a in [0.05, 0.3, 1.0, 2.0, 4.0]:
        x = np.array([0.0, a]); c = np.array([np.cosh(a / 2), -1.0])
        A, PmS, Pm, D, K, md = build_mats(x, a)
        out["primeS_vs_prime_maxdiff"] = max(out["primeS_vs_prime_maxdiff"], md)
        archQ = float(c @ A @ c); primeQ = float(c @ PmS @ c)
        polecorr = 2 * BIGPSI_I2 ** 2 * float(c @ np.sinh(x / 2)) ** 2
        P = primeQ + polecorr; Dv = archQ - P; Kv = float(c @ K @ c)
        ray = Dv / float(c @ c)
        print(f"{a:>5} {archQ:>+13.8f} {primeQ:>+13.8f} {polecorr:>+10.6f} {P:>+13.8f} {Dv:>+11.8f} {Kv:>+12.8f} {ray:>+10.7f} {Dv/abs(archQ):>+9.5f}")
        out["n2"].append({"a": a, "archQ": archQ, "prime_form": primeQ, "pole_corr": polecorr,
                          "primePoleQ": P, "D": Dv, "Kpsi_direct": Kv, "rayleigh_D": ray,
                          "rel_margin_D_over_absA": Dv / abs(archQ),
                          "verdict": "SURVIVES (NUMERICAL)" if Dv >= 0 else "KILLED at this instance (NUMERICAL)"})
    # ---- scan n=3..6 ----
    print("\n-- scan: min eigenvalue of D-form projected on ker l_c (endpoints have an exact-0 duplicate-point mode: see notes) --")
    print(f"{'L':>4} {'config':>14} {'n':>2} {'lam_min_D':>14} {'lam_max_abs_D':>14} {'rel_margin':>11}")
    for L in [0.5, 1.0, 2.0, 4.0]:
        for n in range(3, 7):
            for name, x in configs(L, n).items():
                A, PmS, Pm, D, K, md = build_mats(x, L)
                out["primeS_vs_prime_maxdiff"] = max(out["primeS_vs_prime_maxdiff"], md)
                lam, lamabs, cvec = proj_eigmin(D, x)
                lamK, _, _ = proj_eigmin(K, x)
                print(f"{L:>4} {name:>14} {n:>2} {lam:>+14.9f} {lamabs:>14.6f} {lam/lamabs:>+11.3e}   (K-form proj min {lamK:+.9f})")
                out["scan"].append({"L": L, "config": name, "n": n, "lam_min_D": float(lam),
                                    "lam_max_abs_D": float(lamabs), "lam_min_Kform": float(lamK),
                                    "verdict": "SURVIVES (NUMERICAL)" if lam >= -5e-6 else "KILLED at this instance (NUMERICAL)"})
    print(f"\nprimeTermS(matched S) vs primeTerm: max |diff| over all lags used = {out['primeS_vs_prime_maxdiff']:.3e} (primeTermS_eq check)")

    # ---- mpmath confirmation, dps=40: n=2 table + extreme scan cases ----
    print("\n-- mpmath confirmation (dps=40) --")
    import mpmath as mp
    mp.mp.dps = 40
    def triG_mp(y):
        ay = abs(y)
        if ay <= 1: return mp.mpf(2) / 3 - y * y + ay ** 3 / 2
        if ay <= 2: return (2 - ay) ** 3 / 6
        return mp.mpf(0)
    def prime_mp(d):
        d = mp.mpf(d); nmax = int(mp.e ** (2 + abs(d))) + 2; s = mp.mpf(0)
        for n in range(2, nmax + 1):
            p, _ = _factor_prime(n)
            if p is None: continue
            ln = mp.log(n)
            s += mp.log(p) / mp.sqrt(n) * (triG_mp(ln + d) + triG_mp(-ln + d))
        return s
    def arch_mp(d):
        # oscillatory lags d >= 2 need period subdivision: plain quad at d=4 differs from
        # the subdivided value by 2.1e-6 in archTerm(4) (checked standalone, dps=40)
        d = mp.mpf(d)
        f = lambda t: (mp.mpf(1) if t == 0 else (mp.sin(t / 2) / (t / 2)) ** 2) ** 2 \
            * mp.cos(t * d) * (mp.re(mp.digamma(mp.mpc(mp.mpf("0.25"), t / 2))) - mp.log(mp.pi))
        if d >= 2:
            nper = int(100 * d / mp.pi) + 1
            pts = [mp.mpf(0)] + [mp.mpf(k) * mp.pi / d for k in range(1, nper + 1)] + [mp.inf]
            return (1 / mp.pi) * mp.quad(f, pts)
        return (1 / mp.pi) * mp.quad(f, [0, mp.inf])
    _amp, _pmp = {}, {}
    def Amp(d):
        k = str(d)
        if k not in _amp: _amp[k] = arch_mp(k)
        return _amp[k]
    def Pmp(d):
        k = str(d)
        if k not in _pmp: _pmp[k] = prime_mp(k)
        return _pmp[k]
    bI2 = 16 * mp.sinh(mp.mpf("0.25")) ** 2
    conf = {"n2_mp": []}
    for a in ["0.05", "0.3", "1", "2", "4"]:
        av = mp.mpf(a); A0, Aa = Amp("0"), Amp(a)
        c0, c1 = mp.cosh(av / 2), mp.mpf(-1)
        archQ = c0 * c0 * A0 + 2 * c0 * c1 * Aa + A0
        primeQ = c0 * c0 * Pmp("0") + 2 * c0 * c1 * Pmp(a) + Pmp("0")
        polecorr = 2 * bI2 ** 2 * (c1 * mp.sinh(av / 2)) ** 2
        D = archQ - (primeQ + polecorr)
        print(f"a={a}: archQ={mp.nstr(archQ, 18)} primeForm={mp.nstr(primeQ, 18)} poleCorr={mp.nstr(polecorr, 12)} D={mp.nstr(D, 18)}")
        conf["n2_mp"].append({"a": a, "archQ": mp.nstr(archQ, 20), "prime_form": mp.nstr(primeQ, 20),
                              "pole_corr": mp.nstr(polecorr, 20), "D": mp.nstr(D, 20)})
    def proj_eigmin_mp(xlist):
        n = len(xlist); xs = [mp.mpf(str(v)) for v in xlist]
        ell = [mp.cosh(v / 2) for v in xs]; sv = [mp.sinh(v / 2) for v in xs]
        Dm = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                d = str(abs(float(xs[i] - xs[j])))
                Dm[i, j] = Amp(d) - Pmp(d) - 2 * bI2 ** 2 * sv[i] * sv[j]
        # ker basis: b_k = e_k - (ell_k/ell_0) e_0, then Gram-Schmidt orthonormalisation
        cols = []
        for k in range(1, n):
            v = mp.matrix(n, 1); v[k] = 1; v[0] = -ell[k] / ell[0]; cols.append(v)
        Z = []
        for v in cols:
            w = mp.matrix(v)
            for q in Z:
                w -= q * (q.T * w)[0]
            Z.append(w / mp.sqrt((w.T * w)[0]))
        m = n - 1; Mp_ = mp.matrix(m, m)
        for i in range(m):
            for j in range(m):
                Mp_[i, j] = (Z[i].T * Dm * Z[j])[0]
        ev = mp.eigsy(Mp_, eigvals_only=True)
        return [mp.nstr(v, 18) for v in ev], mp.nstr(min(ev), 18)
    extreme = [("clustered0.1L L=0.5 n=6", list(configs(0.5, 6)["clustered0.1L"])),
               ("clustered0.1L L=0.5 n=5", list(configs(0.5, 5)["clustered0.1L"])),
               ("equispaced L=4 n=6", list(configs(4.0, 6)["equispaced"])),
               ("endpoints L=2 n=4", list(configs(2.0, 4)["endpoints"]))]
    conf["proj_mp"] = []
    for name, x in extreme:
        ev, mn = proj_eigmin_mp(x)
        print(f"{name}: projected D eigenvalues (mp dps40) min = {mn}; all = {ev}")
        conf["proj_mp"].append({"case": name, "lam_min_mp": mn, "eigvals_mp": ev})
    out["mpmath_confirmation"] = conf
    print("\nJSON head:", json.dumps(out)[:200], "...")
    return out

if __name__ == "__main__":
    main()
