#!/usr/bin/env python3
"""
NUMERICAL adversarial kill arm for the open Prop `SemilocalArchDomAll`
(independent arm, Turbo's bench, 2026-11-13). Labels: NUMERICAL only.

=== STATEMENT (source of truth: GPPVerify WeilGramContraction.lean /
WeilGramArchPD.lean / WeilGramSemilocal.lean @ 5fb6c8d0, as transcribed in
research/flint_archdom_scout_2026_11_12.py; this script implements the kernels
from those definitions itself and does NOT import the scout's code) ===

  SemilocalArchDomAll := forall L > 0, SemilocalArchDom (primesBelowExp L) L
  SemilocalArchDom S L := forall n (x c : Fin n -> R),
      (forall i j, |x i - x j| <= L) -> (sum i, c i cosh(x i / 2) = 0) ->
        primePoleQ S x c <= archQ x c
  archQ x c      = kform archTerm x c
  primePoleQ S x c = kform (primeTermS S) x c
                     + 2 * bigPsiI2^2 * (sum i, c i sinh(x i / 2))^2
  primeTermS S d = sum over n with minFac(n) in S of
                     vonMangoldt(n)/sqrt(n) * (triG(log n + d) + triG(-log n + d))
  triG(y) = 2/3 - y^2 + |y|^3/2  (|y| <= 1);  (2-|y|)^3/6  (1 < |y| <= 2);  0 else
  bigPsiI2 = 16 * sinh(1/4)^2
  archTerm(d) = (1/pi) * int_0^inf  [(sin(t/2)/(t/2))^2]^2 cos(t d)
                * (Re digamma(1/4 + i t/2) - log pi) dt

By the PROVED primeTermS_eq, for |d| <= L and S = primesBelowExp L (which
covers L), primeTermS S d = primeTerm d, the unrestricted sum; this script
computes primeTerm. The defect is D(c) := archQ - primePoleQ >= 0 on
ker l_c, l_c(x)_i = cosh(x_i/2). For fixed x the worst c is the projected
minimum eigenvalue of the D-matrix
  D_ij = archTerm(x_i-x_j) - primeTerm(x_i-x_j)
         - 2 bigPsiI2^2 sinh(x_i/2) sinh(x_j/2)
on ker l_c. KILL = one verified negative value (Kill Criterion 1: kills
SemilocalArchDomAll and, by the banked HPsi <-> SemilocalArchDomAll, H_psi).

=== METHOD (fresh implementation, three independent routes) ===
  Route 1 (exact float64): primeTerm by direct prime-power enumeration;
    archTerm by adaptive quadrature with node subdivision + analytic tail.
  Route 2 (surrogate, for global optimisation only): primeTerm from a
    segmented-sieve histogram in u = log n (bin width 1e-4) convolved with
    triG on a d-grid; archTerm from a trapezoid cosine transform; cubic
    splines. Surrogate error is MEASURED against Route 1, never assumed.
  Route 3 (mpmath dps40): independent high-precision kernels for validation
    against Flint's banked values and for final candidate verification.
Assault: scipy differential evolution + Nelder-Mead basin starts (prime-log
clusters, near-coalescent, two-cluster, endpoints), then a global-shift
refinement of the best shapes (the statement constrains only the spread,
so translations are legal and are searched), then exact re-verification
of the smallest candidates with explicit error budgets.
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, "flint_adversarial_kill_2026_11_13_results.json")
CACHE = "/tmp/flint_kill_grids.npz"
SEED = 20261113

LOGPI = math.log(math.pi)
BPSI2 = 16.0 * math.sinh(0.25) ** 2          # bigPsiI2
POLE = 2.0 * BPSI2 ** 2                      # pole coefficient 2*bigPsiI2^2
BANKED_ARCH0 = -1.54293305945377535652       # Flint banked archTerm(0), dps40
BANKED_N2 = {5.0: 0.0132705037417421082, 6.0: 0.0366453810155422438,
             8.0: 0.240200637349104986, 10.0: 1.90613027971460356}  # dps40 D

OUT = {"label": "NUMERICAL", "experiment": "flint_adversarial_kill",
       "statement": "SemilocalArchDomAll (GPPVerify @5fb6c8d0): D(c) = archQ - primePoleQ >= 0 on ker l_c"}


def save_out():
    with open(OUT_JSON, "w") as fh:
        json.dump(OUT, fh, indent=1)


# ================= triG =================
def triG(y):
    ay = abs(y)
    if ay <= 1.0:
        return 2.0 / 3.0 - y * y + ay ** 3 / 2.0
    if ay <= 2.0:
        return (2.0 - ay) ** 3 / 6.0
    return 0.0


def triG_np(y):
    y = np.asarray(y, float)
    ay = np.abs(y)
    out = np.zeros_like(ay)
    m1 = ay <= 1.0
    out[m1] = 2.0 / 3.0 - y[m1] ** 2 + ay[m1] ** 3 / 2.0
    m2 = (~m1) & (ay <= 2.0)
    out[m2] = (2.0 - ay[m2]) ** 3 / 6.0
    return out


# ================= exact float64 prime kernel =================
def _sieve_primes(nmax):
    sieve = np.ones(nmax + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(nmax ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return np.flatnonzero(sieve)


_PP_CACHE = {}


def prime_power_table(nmax):
    """Arrays (u=log n, w=log p / sqrt(n)) over all prime powers n=p^k<=nmax."""
    if nmax in _PP_CACHE:
        return _PP_CACHE[nmax]
    primes = _sieve_primes(nmax)
    us, ws = [], []
    for p in primes:
        lp = math.log(p)
        q, k = int(p), 1
        while q <= nmax:
            us.append(k * lp)
            ws.append(lp / math.sqrt(q))
            if q > nmax // p:
                break
            q *= p
            k += 1
    tab = (np.array(us), np.array(ws))
    _PP_CACHE[nmax] = tab
    return tab


def prime_term_exact(d):
    """Unrestricted primeTerm(d), direct enumeration. Feasible for |d| <= 13
    (table to e^{2+|d|}); raises beyond that (use prime_term_far)."""
    d = abs(float(d))
    nmax = int(math.exp(2.0 + d)) + 3
    if nmax > int(math.exp(15.05)):
        raise ValueError("prime_term_exact: |d| too large, use prime_term_far")
    u, w = prime_power_table(nmax)
    return float(np.sum(w * (triG_np(u + d) + triG_np(u - d))))


# ================= exact float64 arch kernel (adaptive quad) =================
from scipy.integrate import quad as _quad
from scipy.special import digamma as _digamma


def _arch_integrand(t, d):
    if t == 0.0:
        return float(_digamma(0.25).real) - LOGPI
    s = math.sin(t / 2.0) / (t / 2.0)
    return s ** 4 * math.cos(t * d) * (float(_digamma(0.25 + 0.5j * t).real) - LOGPI)


def arch_quad(d, T=300.0):
    """archTerm(|d|) by adaptive quadrature: breakpoints at the integrand's
    own features (sinc zeros 2*pi*k and cosine half-periods pi/d), plus an
    analytic tail bound. Returns (value, error_estimate)."""
    d = abs(float(d))
    pts = set()
    k = 1
    while 2 * math.pi * k < 60.0:
        pts.add(2 * math.pi * k)
        k += 1
    if d > 0.25:
        k = 1
        while k * math.pi / d < 80.0:
            pts.add(k * math.pi / d)
            k += 1
    for v in (1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 80.0, 160.0):
        pts.add(v)
    val, err = _quad(_arch_integrand, 0.0, T, args=(d,), limit=600,
                     epsabs=1e-12, epsrel=1e-12, points=sorted(pts))
    tail = (16.0 / math.pi) * (math.log(T) / (3.0 * T ** 3) + 1.0 / (9.0 * T ** 3))
    return val / math.pi, err / math.pi + tail


# ================= defect matrix / projected eigenvalue =================
def D_matrix(x, archf, primef):
    x = np.asarray(x, float)
    dx = np.abs(x[:, None] - x[None, :])
    A = archf(dx)
    P = primef(dx)
    s = np.sinh(x / 2.0)
    return A - P - POLE * np.outer(s, s)


def proj_lam(D, x):
    """(lam_min, lam_max_abs, worst c) of D restricted to ker l_c, c unit."""
    x = np.asarray(x, float)
    ell = np.cosh(x / 2.0).reshape(1, -1)
    _, _, Vt = np.linalg.svd(ell)
    Z = Vt[1:].T
    M = Z.T @ D @ Z
    w, v = np.linalg.eigh(M)
    return float(w[0]), float(np.max(np.abs(w))), Z @ v[:, 0]


def min_pair_sep(x):
    x = np.sort(np.asarray(x, float))
    if len(x) < 2:
        return float("inf")
    return float(np.min(np.diff(x)))


# ================= far prime kernel (segmented, exact float64) =================
import numba


@numba.njit(cache=True)
def _seg_prime_sum(lo, hi, small_primes, d):
    """Sum of log(p)/sqrt(p) * (triG(log p - d) + triG(log p + d)) over primes
    p in [lo, hi), plus prime powers p^k (k>=2) in range handled by caller.
    Odd-only segmented sieve; lo must be odd."""
    size = (hi - lo + 1) // 2
    is_comp = np.zeros(size, dtype=np.bool_)
    for p in small_primes:
        if p == 2:
            continue
        start = (lo + p - 1) // p * p
        if start < p * p:
            start = p * p
        if start % 2 == 0:
            start += p
        for m in range(start, hi, 2 * p):
            is_comp[(m - lo) // 2] = True
    total = 0.0
    for i in range(size):
        if not is_comp[i]:
            n = lo + 2 * i
            if n >= 3:
                u = math.log(n)
                t = 0.0
                a1 = abs(u - d)
                if a1 <= 1.0:
                    t += 2.0 / 3.0 - (u - d) ** 2 + a1 ** 3 / 2.0
                elif a1 <= 2.0:
                    t += (2.0 - a1) ** 3 / 6.0
                a2 = abs(u + d)
                if a2 <= 1.0:
                    t += 2.0 / 3.0 - (u + d) ** 2 + a2 ** 3 / 2.0
                elif a2 <= 2.0:
                    t += (2.0 - a2) ** 3 / 6.0
                total += (u / math.sqrt(n)) * t
    return total


def prime_term_far(d):
    """Exact primeTerm(d) for any d >= 0 by segmented sieve over the only
    window that contributes, n in [e^{d-2}, e^{d+2}] (triG support), plus all
    prime powers. Pairwise (segmented) accumulation; float64."""
    d = abs(float(d))
    lo = max(3, int(math.exp(d - 2.0)) - 2)
    if lo % 2 == 0:
        lo += 1
    hi = int(math.exp(d + 2.0)) + 4
    small = _sieve_primes(int(hi ** 0.5) + 2)
    seg = 1 << 25  # even, so every segment start stays odd
    parts = []
    for a in range(lo, hi, seg):
        b = min(hi, a + seg)
        parts.append(_seg_prime_sum(a, b, small.astype(np.int64), d))
    # prime powers p^k, k >= 2, n <= e^{d+2} (primes themselves already counted
    # only inside [lo, hi); small primes below lo contribute via triG(u+d)
    # when d < 2, handled by extending the window: for d < 2, lo is ~3 anyway)
    total = math.fsum(parts)
    nmax = int(math.exp(d + 2.0)) + 3
    for p in _sieve_primes(int(nmax ** 0.5) + 2):
        lp = math.log(int(p))
        q = int(p) * int(p)
        while q <= nmax:
            u = math.log(q)
            total += (lp / math.sqrt(q)) * (triG(u - d) + triG(u + d))
            if q > nmax // p:
                break
            q *= p
    if d < 2.0:
        # small primes below the segmented window (n < lo) with triG(u+d) mass
        for n, p in ((2, 2),):
            u = math.log(n)
            total += (math.log(p) / math.sqrt(n)) * (triG(u - d) + triG(u + d))
    return float(total)


# ================= mpmath dps40 kernels =================
def mp_kernels(dps=40):
    import mpmath as mp
    mp.mp.dps = dps

    def triG_mp(y):
        ay = abs(y)
        if ay <= 1:
            return mp.mpf(2) / 3 - y * y + ay ** 3 / 2
        if ay <= 2:
            return (2 - ay) ** 3 / 6
        return mp.mpf(0)

    _pp = {}

    def pp_table(nmax):
        if nmax not in _pp:
            sieve = bytearray(b"\x01") * (nmax + 1)
            sieve[0:2] = b"\x00\x00"
            for i in range(2, int(nmax ** 0.5) + 1):
                if sieve[i]:
                    sieve[i * i::i] = b"\x00" * len(sieve[i * i::i])
            tab = []
            for p in range(2, nmax + 1):
                if sieve[p]:
                    q = p
                    while q <= nmax:
                        tab.append((q, p))
                        if q > nmax // p:
                            break
                        q *= p
            _pp[nmax] = tab
        return _pp[nmax]

    def prime_mp(d):
        d = mp.mpf(d)
        nmax = int(mp.e ** (2 + abs(d))) + 2
        s = mp.mpf(0)
        for n, p in pp_table(nmax):
            ln = mp.log(n)
            s += mp.log(p) / mp.sqrt(n) * (triG_mp(ln + d) + triG_mp(ln - d))
        return s

    def arch_mp(d):
        d = mp.mpf(d)
        f = lambda t: (mp.mpf(1) if t == 0 else (mp.sin(t / 2) / (t / 2)) ** 2) ** 2 \
            * mp.cos(t * d) * (mp.re(mp.digamma(mp.mpc(mp.mpf("0.25"), t / 2))) - mp.log(mp.pi))
        if d >= 2:
            nper = int(100 * d / mp.pi) + 1
            pts = [mp.mpf(0)] + [mp.mpf(k) * mp.pi / d for k in range(1, nper + 1)] \
                + [mp.inf]
            return (1 / mp.pi) * mp.quad(f, pts)
        return (1 / mp.pi) * mp.quad(f, [0, mp.inf])

    return mp, prime_mp, arch_mp


def stage_validate():
    """Route 3 validation against Flint's banked dps40 values + Route 1
    cross-checks (arch quad vs banked, prime exact vs segmented far)."""
    t0 = time.time()
    mp, prime_mp, arch_mp = mp_kernels(40)
    res = {"label": "NUMERICAL"}
    a0 = arch_mp("0")
    res["archTerm0_mp"] = mp.nstr(a0, 25)
    res["archTerm0_banked"] = repr(BANKED_ARCH0)
    res["archTerm0_absdiff"] = float(abs(a0 - mp.mpf(repr(BANKED_ARCH0))))
    aq, ae = arch_quad(0.0)
    res["archTerm0_quad"] = aq
    res["archTerm0_quad_err_est"] = ae
    res["archTerm0_quad_vs_mp_absdiff"] = float(abs(mp.mpf(aq) - a0))
    bI2 = 16 * mp.sinh(mp.mpf("0.25")) ** 2
    rows = []
    cache = {}

    def Amp(s):
        if s not in cache:
            cache[s] = arch_mp(s)
        return cache[s]

    def Pmp(s):
        k = "p" + s
        if k not in cache:
            cache[k] = prime_mp(s)
        return cache[k]

    for a, banked in sorted(BANKED_N2.items()):
        A0, Aa = Amp("0"), Amp(repr(a))
        P0, Pa = Pmp("0"), Pmp(repr(a))
        av = mp.mpf(repr(a))
        c0, c1 = mp.cosh(av / 2), mp.mpf(-1)
        archQ = c0 * c0 * A0 + 2 * c0 * c1 * Aa + A0
        primeQ = c0 * c0 * P0 + 2 * c0 * c1 * Pa + P0
        pole = 2 * bI2 ** 2 * (c1 * mp.sinh(av / 2)) ** 2
        D = archQ - (primeQ + pole)
        # second route: K-form matrix (identity on ker l_c, PROVED)
        K00 = 2 * bI2 ** 2 - P0 + A0
        Kaa = 2 * bI2 ** 2 * mp.cosh(av / 2) - Pa + Aa
        DK = c0 * c0 * K00 + 2 * c0 * c1 * Kaa + K00
        rows.append({"a": a, "D_mp": mp.nstr(D, 20), "D_banked": repr(banked),
                     "absdiff_vs_banked": float(abs(D - mp.mpf(repr(banked)))),
                     "D_Kform_mp": mp.nstr(DK, 20),
                     "Kform_vs_D_absdiff": float(abs(D - DK))})
        print(f"validate n=2 a={a}: D={mp.nstr(D, 15)} banked={banked:.15f} "
              f"diff={rows[-1]['absdiff_vs_banked']:.3e} Kform_diff={rows[-1]['Kform_vs_D_absdiff']:.3e}")
    res["n2_rows"] = rows
    # Route1 vs Route-far prime cross-check at moderate lags
    cc = []
    for d in (3.0, 7.0, 11.0):
        e1 = prime_term_exact(d)
        e2 = prime_term_far(d)
        cc.append({"d": d, "exact": e1, "far": e2, "absdiff": abs(e1 - e2)})
        print(f"prime cross-check d={d}: exact={e1:.12f} far={e2:.12f} diff={abs(e1-e2):.3e}")
    res["prime_crosschecks"] = cc
    res["runtime_s"] = time.time() - t0
    OUT["validation"] = res
    save_out()
    print(f"validation done in {res['runtime_s']:.1f}s; arch0 diff vs banked = {res['archTerm0_absdiff']:.3e}")
    return res


# ================= surrogate grids (Route 2) =================
DGRID_MAX = 20.0
DGRID_STEP = 0.005
DU = 1e-4
U_MAX = 22.02
NMAX_HIST = int(math.exp(U_MAX)) + 10

_splines = {}


@numba.njit(cache=True)
def _hist_segment(lo, hi, small_primes, Hl, DU):
    size = (hi - lo + 1) // 2
    is_comp = np.zeros(size, dtype=np.bool_)
    for p in small_primes:
        if p == 2:
            continue
        start = (lo + p - 1) // p * p
        if start < p * p:
            start = p * p
        if start % 2 == 0:
            start += p
        for m in range(start, hi, 2 * p):
            is_comp[(m - lo) // 2] = True
    for i in range(size):
        if not is_comp[i]:
            n = lo + 2 * i
            if n >= 3:
                u = math.log(float(n))
                Hl[int(u / DU)] += u / math.sqrt(float(n))


def build_histogram():
    nbins = int(U_MAX / DU) + 2
    H = np.zeros(nbins)
    small = _sieve_primes(int(NMAX_HIST ** 0.5) + 2).astype(np.int64)
    seg = 1 << 26
    lo = 3
    while lo <= NMAX_HIST:
        hi = min(NMAX_HIST + 1, lo + seg)
        _hist_segment(lo, hi, small, H, DU)
        lo = hi if hi % 2 == 1 else hi + 1
    # prime 2 and prime powers p^k, k >= 2
    u2 = math.log(2.0)
    H[int(u2 / DU)] += u2 / math.sqrt(2.0)
    for p in _sieve_primes(int(NMAX_HIST ** 0.5) + 2):
        lp = math.log(int(p))
        q = int(p) * int(p)
        while q <= NMAX_HIST:
            u = math.log(q)
            H[int(u / DU)] += lp / math.sqrt(q)
            if q > NMAX_HIST // p:
                break
            q *= p
    return H


def prime_from_hist(H, du, dgrid):
    centers = (np.arange(len(H)) + 0.5) * du
    out = np.empty(len(dgrid))
    for i, d in enumerate(dgrid):
        m = (centers >= d - 2.0) & (centers <= d + 2.0)
        F = float(np.dot(H[m], triG_np(centers[m] - d)))
        m2 = centers <= 2.0 - d
        Fm = float(np.dot(H[m2], triG_np(centers[m2] + d))) if d < 2.0 else 0.0
        out[i] = F + Fm
    return out


def build_arch_grid(dgrid):
    from scipy.special import digamma as dg
    t = np.arange(0.0, 300.0 + 1e-9, 0.002)
    s = np.ones_like(t)
    nz = t > 0
    s[nz] = np.sin(t[nz] / 2.0) / (t[nz] / 2.0)
    g = s ** 4 * (dg(0.25 + 0.5j * t).real - LOGPI)
    g[0] = float(dg(0.25).real) - LOGPI
    w = np.full_like(t, 0.002)
    w[0] = w[-1] = 0.001
    gw = g * w
    vals = np.empty(len(dgrid))
    for a in range(0, len(dgrid), 96):
        C = np.cos(np.outer(dgrid[a:a + 96], t))
        vals[a:a + 96] = C @ gw / math.pi
    return vals


def ensure_grids():
    from scipy.interpolate import CubicSpline
    if _splines:
        return _splines
    if os.path.exists(CACHE):
        z = np.load(CACHE)
        dgrid, pg, ag = z["dgrid"], z["prime_grid"], z["arch_grid"]
        meta = {"source": "cache"}
    else:
        t0 = time.time()
        H = build_histogram()
        hist_s = time.time() - t0
        dgrid = np.arange(0.0, DGRID_MAX + 1e-9, DGRID_STEP)
        pg = prime_from_hist(H, DU, dgrid)
        ag = build_arch_grid(dgrid)
        np.savez(CACHE, dgrid=dgrid, prime_grid=pg, arch_grid=ag, H=H)
        meta = {"source": "built", "histogram_s": hist_s}
    _splines["dgrid"] = dgrid
    _splines["prime"] = CubicSpline(dgrid, pg)
    _splines["arch"] = CubicSpline(dgrid, ag)
    _splines["meta"] = meta
    return _splines


def arch_sur(dx):
    return _splines["arch"](np.abs(dx))


def prime_sur(dx):
    return _splines["prime"](np.abs(dx))


def stage_grids():
    t0 = time.time()
    sp = ensure_grids()
    rng = np.random.default_rng(SEED)
    res = {"label": "NUMERICAL", "meta": sp["meta"]}
    # prime surrogate error vs exact (direct enumeration, d <= 12)
    errs = []
    for d in rng.uniform(0.5, 12.0, 18):
        errs.append(abs(float(prime_sur(d)) - prime_term_exact(float(d))))
    res["prime_spline_vs_exact_maxerr_d_le_12"] = float(np.max(errs))
    # prime surrogate error vs segmented far route (independent path)
    errs_far = []
    for d in (13.5, 16.0, 18.5, 19.8):
        errs_far.append(abs(float(prime_sur(d)) - prime_term_far(d)))
    res["prime_spline_vs_far_absdiffs"] = [float(e) for e in errs_far]
    # histogram binning error: DU vs 2*DU at sample grid points
    z = np.load(CACHE)
    H = z["H"]
    H2 = H[:len(H) // 2 * 2].reshape(-1, 2).sum(axis=1)
    sample = np.array([2.0, 5.0, 9.0, 13.0, 17.0, 19.5])
    p1 = prime_from_hist(H, DU, sample)
    p2 = prime_from_hist(H2, 2 * DU, sample)
    res["hist_binning_diff_DU_vs_2DU"] = [float(abs(a - b)) for a, b in zip(p1, p2)]
    # arch surrogate error vs adaptive quad
    aerrs = []
    for d in rng.uniform(0.0, 20.0, 14):
        v, _ = arch_quad(float(d))
        aerrs.append(abs(float(arch_sur(d)) - v))
    res["arch_spline_vs_quad_maxerr"] = float(np.max(aerrs))
    # lambda-level surrogate error on random configs
    lerrs = []
    for (n, L) in ((3, 2.0), (4, 5.0), (5, 10.0), (6, 10.0)):
        x = np.sort(rng.uniform(-L / 2, L / 2, n))
        lam_s, _, _ = proj_lam(D_matrix(x, arch_sur, prime_sur), x)
        lam_e, _, _ = proj_lam(D_matrix(x, lambda dx: np.vectorize(lambda v: arch_quad(v)[0])(dx),
                                        lambda dx: np.vectorize(prime_term_exact)(dx)), x)
        lerrs.append({"n": n, "L": L, "lam_surrogate": lam_s, "lam_exact": lam_e,
                      "absdiff": abs(lam_s - lam_e)})
    res["lambda_level_checks"] = lerrs
    res["runtime_s"] = time.time() - t0
    OUT["grids"] = res
    save_out()
    print("grids:", json.dumps({k: v for k, v in res.items() if k != 'lambda_level_checks'}, default=str)[:600])
    print("lambda checks:", lerrs)
    return res


# ================= assault (Route 2 search) =================
from scipy.optimize import differential_evolution, minimize, minimize_scalar

EVALS = {"count": 0}
CANDIDATES = []


def lam_sur(x):
    EVALS["count"] += 1
    x = np.asarray(x, float)
    D = D_matrix(x, arch_sur, prime_sur)
    lam, _, _ = proj_lam(D, x)
    return lam


def lam_sur_penalised(x, L):
    x = np.asarray(x, float)
    spread = float(np.max(x) - np.min(x)) if len(x) else 0.0
    pen = 0.0
    if spread > L:
        pen += 1e4 * (spread - L) ** 2 + 10.0 * (spread - L)
    if np.max(np.abs(x)) > 70.0:
        pen += 1e3 * (np.max(np.abs(x)) - 70.0)
    return lam_sur(x) + pen


def adversarial_starts(n, L, rng):
    starts = []
    half = (n - 1) / 2.0
    starts.append(("equispaced", np.linspace(-L / 2, L / 2, n)))
    starts.append(("endpoints", np.array([-L / 2 if i % 2 == 0 else L / 2 for i in range(n)], float)))
    for eps in (1e-3, 3e-2, 0.3):
        starts.append((f"coalescent_eps{eps}", (np.arange(n) - half) * eps))
    for k in (n // 2,):
        for a in (L / 4, 0.45 * L):
            x = np.array([-a] * k + [a] * (n - k), float)
            x = x + rng.normal(0, 1e-3, n)
            starts.append((f"twocluster_a{a:.2f}", np.sort(x)))
    primes = _sieve_primes(10000)
    lp = np.log(primes.astype(float))
    for k0 in (0, 25, 100, 400, 1200):
        if k0 + n <= len(lp):
            raw = lp[k0:k0 + n]
            starts.append((f"logprime_raw_k{k0}", raw.copy()))
            starts.append((f"logprime_centred_k{k0}", raw - np.mean(raw)))
            starts.append((f"logprime_leftedge_k{k0}", raw - raw[0] - L / 2))
    for r in range(3):
        starts.append((f"random{r}", np.sort(rng.uniform(-L / 2, L / 2, n))))
    return starts


def run_assault():
    ensure_grids()
    rng = np.random.default_rng(SEED)
    per = []
    for L in (2.0, 5.0, 10.0, 20.0):
        for n in range(3, 11):
            t0 = time.time()
            evals0 = EVALS["count"]
            best = {"lam": float("inf")}
            best_nd = {"lam": float("inf")}

            def consider(x, source):
                x = np.asarray(x, float)
                lam = lam_sur(x)
                rec = {"lam": float(lam), "x": [float(v) for v in x],
                       "sep": min_pair_sep(x), "source": source, "n": n, "L": L}
                if lam < best["lam"]:
                    best.update(rec)
                if rec["sep"] >= 0.02 and lam < best_nd["lam"]:
                    best_nd.update(rec)
                return lam

            maxiter = 300 if n <= 5 else (250 if n <= 8 else 200)
            de = differential_evolution(
                lambda v: lam_sur(v), [(-L / 2, L / 2)] * n,
                seed=SEED + n * 100 + int(L), popsize=10, maxiter=maxiter,
                mutation=(0.3, 1.2), recombination=0.8, polish=False, tol=0)
            consider(de.x, "DE")
            nm0 = minimize(lambda v: lam_sur_penalised(v, L), de.x,
                           method="Nelder-Mead",
                           options={"maxiter": 2500, "xatol": 1e-4, "fatol": 1e-10})
            consider(nm0.x, "DE+NM")
            for name, x0 in adversarial_starts(n, L, rng):
                consider(x0, f"start:{name}")
                r = minimize(lambda v: lam_sur_penalised(v, L), x0,
                             method="Nelder-Mead",
                             options={"maxiter": 2500, "xatol": 1e-4, "fatol": 1e-10})
                consider(r.x, f"NM:{name}")
            row = {"n": n, "L": L, "best": best, "best_nondegenerate": best_nd,
                   "evals": EVALS["count"] - evals0, "runtime_s": time.time() - t0}
            per.append(row)
            CANDIDATES.append(best)
            if best_nd["lam"] < float("inf"):
                CANDIDATES.append(best_nd)
            print(f"L={L} n={n}: best lam={best['lam']:+.6e} (sep {best['sep']:.4f}, {best['source']}); "
                  f"nondeg lam={best_nd['lam']:+.6e} (sep {best_nd['sep']:.4f}, {best_nd['source']}); "
                  f"evals={row['evals']}")
    # ---- global-shift refinement of the best shapes ----
    shifted = []
    seen = set()
    for cand in sorted(CANDIDATES, key=lambda c: c["lam"])[:24]:
        key = tuple(np.round(np.sort(cand["x"]) - np.mean(cand["x"]), 2))
        if key in seen:
            continue
        seen.add(key)
        x0 = np.asarray(cand["x"], float)
        shape = x0 - np.mean(x0)
        ss = np.linspace(-15.0, 15.0, 121)
        vals = [lam_sur(shape + s) for s in ss]
        s_best = float(ss[int(np.argmin(vals))])
        r = minimize_scalar(lambda s: lam_sur(shape + s), bounds=(s_best - 1.0, s_best + 1.0),
                            method="bounded", options={"xatol": 1e-3})
        lam_r = float(r.fun)
        base_lam = float(lam_sur(x0))
        rec = {"lam": min(lam_r, base_lam),
               "x": [float(v) for v in (shape + (r.x if lam_r <= base_lam else np.mean(x0)))],
               "sep": min_pair_sep(shape), "source": "shift_refined",
               "n": cand["n"], "L": cand["L"], "shift": float(r.x if lam_r <= base_lam else np.mean(x0))}
        shifted.append(rec)
        print(f"shift refine n={cand['n']} L={cand['L']}: base {base_lam:+.6e} -> {rec['lam']:+.6e} at shift {rec['shift']:+.2f}")
    allc = CANDIDATES + shifted
    allc.sort(key=lambda c: c["lam"])
    OUT["assault"] = {"label": "NUMERICAL", "per_nL": per,
                      "total_evals": EVALS["count"],
                      "candidates_sorted": allc[:20]}
    save_out()
    return allc


# ================= exact verification =================
def exact_lam(x):
    """Projected lam_min with exact float64 kernels (prime: enumeration for
    lags <= 13, segmented far route above; arch: adaptive quad) and a
    Gershgorin error budget from the quadrature error estimates."""
    x = np.asarray(x, float)
    n = len(x)
    dx = np.abs(x[:, None] - x[None, :])
    A = np.zeros((n, n))
    E = np.zeros((n, n))
    P = np.zeros((n, n))
    acache, pcache = {}, {}
    for i in range(n):
        for j in range(n):
            d = float(dx[i, j])
            if d not in acache:
                acache[d] = arch_quad(d)
            if d not in pcache:
                pcache[d] = prime_term_exact(d) if d <= 13.0 else prime_term_far(d)
            A[i, j], E[i, j] = acache[d]
            P[i, j] = pcache[d]
    s = np.sinh(x / 2.0)
    D = A - P - POLE * np.outer(s, s)
    lam, lamabs, c = proj_lam(D, x)
    err = float(np.max(np.sum(np.abs(E), axis=1))) + 1e-6  # + prime float slack
    return {"lam_exact": float(lam), "err_budget": err, "c_worst": [float(v) for v in c],
            "D_of_c": float(c @ D @ c)}


def exact_lam_fast(x):
    """exact_lam, but prime lags evaluated in one segmented pass when the
    largest lag exceeds 13 (per-lag segmented calls would repeat the sieve).
    Includes a spot check of the multi-lag pass against prime_term_far at
    the largest lag."""
    x = np.asarray(x, float)
    n = len(x)
    dx = np.abs(x[:, None] - x[None, :])
    uniq = sorted(set(float(v) for v in dx.ravel()))
    if max(uniq) > 13.0:
        pmap = prime_terms_multi(uniq)
        spot = abs(pmap[max(uniq)] - prime_term_far(max(uniq)))
    else:
        pmap = {d: prime_term_exact(d) for d in uniq}
        spot = 0.0
    A = np.zeros((n, n))
    E = np.zeros((n, n))
    P = np.zeros((n, n))
    acache = {}
    for i in range(n):
        for j in range(n):
            d = float(dx[i, j])
            if d not in acache:
                acache[d] = arch_quad(d)
            A[i, j], E[i, j] = acache[d]
            P[i, j] = pmap[d]
    s = np.sinh(x / 2.0)
    D = A - P - POLE * np.outer(s, s)
    lam, lamabs, c = proj_lam(D, x)
    err = float(np.max(np.sum(np.abs(E), axis=1))) + 1e-6
    return {"lam_exact": float(lam), "err_budget": err, "c_worst": [float(v) for v in c],
            "D_of_c": float(c @ D @ c), "prime_multi_spotcheck_absdiff": float(spot)}


def mp_verify_candidate(x):
    """mpmath dps40 re-assembly for a candidate whose lags are <= 12
    (prime by enumeration at dps40, arch by subdivided mp quadrature)."""
    mp, prime_mp, arch_mp = mp_kernels(40)
    x = [mp.mpf(repr(float(v))) for v in x]
    n = len(x)
    bI2 = 16 * mp.sinh(mp.mpf("0.25")) ** 2
    ac, pc = {}, {}

    def Amp(d):
        k = mp.nstr(d, 15)
        if k not in ac:
            ac[k] = arch_mp(d)
        return ac[k]

    def Pmp(d):
        k = mp.nstr(d, 15)
        if k not in pc:
            pc[k] = prime_mp(d)
        return pc[k]

    Dm = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            d = abs(x[i] - x[j])
            Dm[i, j] = Amp(d) - Pmp(d) - 2 * bI2 ** 2 * mp.sinh(x[i] / 2) * mp.sinh(x[j] / 2)
    ell = [mp.cosh(v / 2) for v in x]
    cols = []
    for k in range(1, n):
        v = mp.matrix(n, 1)
        v[k] = 1
        v[0] = -ell[k] / ell[0]
        cols.append(v)
    Z = []
    for v in cols:
        w = mp.matrix(v)
        for q in Z:
            w -= q * (q.T * w)[0]
        Z.append(w / mp.sqrt((w.T * w)[0]))
    m = n - 1
    M = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            M[i, j] = (Z[i].T * Dm * Z[j])[0]
    ev = mp.eigsy(M, eigvals_only=True)
    return {"lam_min_mp": mp.nstr(min(ev), 20), "eigvals_mp": [mp.nstr(v, 15) for v in ev]}


def stage_verify():
    if "assault" not in OUT:
        run_assault()
    cands = OUT["assault"]["candidates_sorted"]
    picked, seen = [], set()

    def _add(c):
        key = (c["n"], c["L"], tuple(np.round(np.sort(c["x"]), 2)))
        if key not in seen:
            seen.add(key)
            picked.append(c)

    for c in cands[:4]:
        _add(c)
    per = {(r["n"], r["L"]): r for r in OUT["assault"]["per_nL"]}
    for nl in ((9, 2.0), (5, 2.0), (4, 5.0)):
        if nl in per:
            _add(per[nl]["best"])
    rows = []
    killed = None
    for c in picked:
        x = np.array(c["x"], float)
        r = exact_lam_fast(x)
        rec = {"n": c["n"], "L": c["L"], "source": c["source"], "sep": c["sep"],
               "lam_surrogate": c["lam"], "x": c["x"], **r}
        if r["lam_exact"] < -r["err_budget"]:
            rec["status"] = "NEGATIVE BEYOND ERROR BUDGET -- escalating to mpmath"
            if max(abs(x[i] - x[j]) for i in range(len(x)) for j in range(len(x))) <= 12.0:
                rec["mpmath"] = mp_verify_candidate(x)
            killed = rec
        elif r["lam_exact"] < 0:
            rec["status"] = "negative inside error budget (near-miss, NOT a kill)"
        else:
            rec["status"] = "positive under exact kernels"
        rows.append(rec)
        print(f"verify n={c['n']} L={c['L']} ({c['source']}): surrogate {c['lam']:+.6e} -> "
              f"exact {r['lam_exact']:+.6e} +- {r['err_budget']:.2e} : {rec['status']}")
    OUT["verification"] = {"label": "NUMERICAL", "rows": rows,
                           "killed": killed is not None}
    save_out()
    return rows


def stage_n2ext():
    """n=2 family beyond Flint's range: a in {12, 15, 20}, float-exact
    assembly with explicit error budget (arch by adaptive quad; prime by
    enumeration / segmented far route; K-form as a second assembly)."""
    rows = []
    for a in (12.0, 15.0, 20.0):
        A0, e0 = arch_quad(0.0)
        Aa, ea = arch_quad(a)
        P0 = prime_term_exact(0.0)
        Pa = prime_term_exact(a) if a <= 13.0 else prime_term_far(a)
        c0, c1 = math.cosh(a / 2.0), -1.0
        archQ = c0 * c0 * A0 + 2 * c0 * c1 * Aa + A0
        primeQ = c0 * c0 * P0 + 2 * c0 * c1 * Pa + P0
        pole = POLE * (c1 * math.sinh(a / 2.0)) ** 2
        D = archQ - (primeQ + pole)
        K00 = POLE - P0 + A0
        Kaa = POLE * math.cosh(a / 2.0) - Pa + Aa
        DK = c0 * c0 * K00 + 2 * c0 * c1 * Kaa + K00
        err = (c0 * c0 + 1.0) * e0 + 2 * c0 * ea + 1e-6 * (2 * c0 + c0 * c0)
        rows.append({"a": a, "D": float(D), "D_Kform": float(DK),
                     "err_budget": float(err), "rel_margin_D_over_absA": float(D / abs(archQ)),
                     "verdict": "SURVIVES (NUMERICAL)" if D > err else
                                ("KILLED at this instance (NUMERICAL)" if D < -err else
                                 "INDETERMINATE within error budget")})
        print(f"n2ext a={a}: D={D:+.6f} Kform={DK:+.6f} err=+/-{err:.2e} rel={rows[-1]['rel_margin_D_over_absA']:+.3e}")
    OUT["n2_extension_beyond_10"] = {"label": "NUMERICAL", "rows": rows}
    save_out()
    return rows


def stage_core():
    """Bounded replay core (< 300 s): anchors + exact re-evaluation of the
    recorded sub-L=20 verified candidates + n2ext a=12. The L=20 exact
    verifications (multi-lag segmented pass) are reproduced by the verify
    stage, which exceeds the runner's 300 s cap."""
    aq, ae = arch_quad(0.0)
    print(f"core anchor archTerm(0) = {aq:.15f} (banked {BANKED_ARCH0}, quad err est {ae:.1e})")
    print(f"core anchor primeTerm(0) = {prime_term_exact(0.0):.15f}; "
          f"prime far-vs-exact at d=7: {abs(prime_term_far(7.0) - prime_term_exact(7.0)):.2e}")
    per = {(r["n"], r["L"]): r for r in OUT["assault"]["per_nL"]}
    rows = []
    for nl in ((9, 2.0), (5, 2.0), (4, 5.0), (10, 10.0)):
        cand = per[nl]["best"]
        r = exact_lam(np.array(cand["x"], float))
        rows.append({"n": nl[0], "L": nl[1], "lam_exact": r["lam_exact"],
                     "err_budget": r["err_budget"]})
        print(f"core verify n={nl[0]} L={nl[1]}: lam_exact {r['lam_exact']:+.6e} +- {r['err_budget']:.2e}")
    a = 12.0
    A0, e0 = arch_quad(0.0)
    Aa, ea = arch_quad(a)
    P0, Pa = prime_term_exact(0.0), prime_term_exact(a)
    c0 = math.cosh(a / 2.0)
    D = (c0 * c0 * A0 - 2 * c0 * Aa + A0) - (c0 * c0 * P0 - 2 * c0 * Pa + P0) \
        - POLE * math.sinh(a / 2.0) ** 2
    print(f"core n2ext a=12: D = {D:+.6f}")
    OUT["core_replay"] = {"label": "NUMERICAL", "rows": rows, "n2ext_a12_D": float(D)}
    save_out()
    return rows


STAGES = {"validate": stage_validate, "grids": stage_grids,
          "assault": run_assault, "verify": stage_verify, "n2ext": stage_n2ext,
          "core": stage_core}




# ================= multi-lag exact prime evaluation =================
@numba.njit(cache=True)
def _seg_prime_multi(lo, hi, small_primes, lags, acc):
    """Accumulate primeTerm(lags[k]) into acc[k] over primes in [lo, hi)."""
    size = (hi - lo + 1) // 2
    is_comp = np.zeros(size, dtype=np.bool_)
    for p in small_primes:
        if p == 2:
            continue
        start = (lo + p - 1) // p * p
        if start < p * p:
            start = p * p
        if start % 2 == 0:
            start += p
        for m in range(start, hi, 2 * p):
            is_comp[(m - lo) // 2] = True
    nl = len(lags)
    for i in range(size):
        if not is_comp[i]:
            n = lo + 2 * i
            if n >= 3:
                u = math.log(float(n))
                w = u / math.sqrt(float(n))
                # binary search: first lag >= u - 2
                a, b = 0, nl
                while a < b:
                    mid = (a + b) // 2
                    if lags[mid] < u - 2.0:
                        a = mid + 1
                    else:
                        b = mid
                k = a
                while k < nl and lags[k] <= u + 2.0:
                    d = lags[k]
                    t = 0.0
                    a1 = abs(u - d)
                    if a1 <= 1.0:
                        t += 2.0 / 3.0 - (u - d) ** 2 + a1 ** 3 / 2.0
                    elif a1 <= 2.0:
                        t += (2.0 - a1) ** 3 / 6.0
                    a2 = abs(u + d)
                    if a2 <= 1.0:
                        t += 2.0 / 3.0 - (u + d) ** 2 + a2 ** 3 / 2.0
                    elif a2 <= 2.0:
                        t += (2.0 - a2) ** 3 / 6.0
                    acc[k] += w * t
                    k += 1


def prime_terms_multi(lag_list):
    """Exact primeTerm at many lags in ONE segmented pass over the union
    window [3, e^{max lag + 2}] (prime powers added separately). Returns
    dict lag -> value."""
    lags = np.array(sorted(set(float(abs(l)) for l in lag_list)), float)
    acc = np.zeros(len(lags))
    dmax = float(lags[-1]) if len(lags) else 0.0
    hi = int(math.exp(dmax + 2.0)) + 4
    small = _sieve_primes(int(hi ** 0.5) + 2).astype(np.int64)
    seg = 1 << 26
    lo = 3
    while lo < hi:
        b = min(hi, lo + seg)
        _seg_prime_multi(lo, b, small, lags, acc)
        lo = b if b % 2 == 1 else b + 1
    out = {float(l): float(acc[k]) for k, l in enumerate(lags)}
    # prime 2 and prime powers p^k (k >= 2)
    extras = [(math.log(2.0), math.log(2.0) / math.sqrt(2.0))]
    nmax = hi
    for p in _sieve_primes(int(nmax ** 0.5) + 2):
        lp = math.log(int(p))
        q = int(p) * int(p)
        while q <= nmax:
            extras.append((math.log(q), lp / math.sqrt(q)))
            if q > nmax // p:
                break
            q *= p
    for u, w in extras:
        for l in lags:
            out[float(l)] += w * (triG(u - l) + triG(u + l))
    return out


if __name__ == "__main__":
    todo = sys.argv[1:] or ["validate", "grids", "assault", "verify", "n2ext"]
    if os.path.exists(OUT_JSON):
        OUT.update(json.load(open(OUT_JSON)))
    for name in todo:
        print(f"=== STAGE {name} ===", flush=True)
        STAGES[name]()
    print("ALL STAGES DONE")
