#!/usr/bin/env python3
"""
NUMERICAL F_theta control for the semilocal defect apparatus
(Flint / Discovery 3, 2026-11-13; queue item 4(ii), Kill Criterion 2).

Labels: NUMERICAL only. This is a CONTROL experiment: it validates (or
indicts) the search machinery of flint_archdom_scout /
flint_archdom_margin_extension / flint_primelog_cluster, not a Prop.

=== THE F_theta DATUM, PINNED FROM THE RECORD BEFORE ANY NUMBER IS TRUSTED ===
Source: ARP hunt note (Turbo's workspace, smatrix-lane/ARP-HUNT-2026-11-10.md),
banked entries of 2026-11-10 21:45-22:00 EST, Arb-verified at 256-bit there:
  * The tree's F_theta at theta = 1/4 is Ftheta = xi(s-1/4) * xi(s+1/4).
  * `primeSideFtheta` with weights n^{1/4} + n^{-1/4} is the correct
    finite-place pairing for Ftheta under the tree's `primeSide`
    normalisation: per prime power n = p^k, the pairing weight is
        Lambda(n) * n^{-1/2} * (n^{1/4} + n^{-1/4})
    in place of zeta's Lambda(n) * n^{-1/2}
    (n^{-1/4} = n^{-1/2} * n^{1/4} from xi(s-1/4); n^{-3/4} = n^{-1/2} *
    n^{-1/4} from xi(s+1/4)). Banked constant: at n = 2 the weight factor is
    2^{1/4} + 2^{-1/4} = 2.0301035302564356... > 2 (a_p > 2: the E4 refusal).
  * The same record rules: the prime-side distinction alone is sufficient
    to discriminate F_theta for the hard tests (BT3 / #561); the archimedean
    difference W_inf^{F_theta} != W_inf^{zeta} is real but is a separate
    named gap at the mu/ARP level and must not gate the hard tests.
  * F_theta's banked failures are prime-sector: lambda_min(P-) at K=2
    = -6.9440067888e-5; E4 refusal at the (p, k=1) coefficients.

TRANSPORT INTO THE SEMILOCAL APPARATUS. The semilocal prime kernel of
flint_archdom_scout is primeTerm(d) = primeSide(h_d) for the even profile
h_d(u) = triG(u + d) + triG(-u + d):
    primeTerm(d) = sum_{n = p^k} (log p) n^{-1/2} (triG(log n + d) + triG(-log n + d)).
Substituting the pinned F_theta pairing weights gives the control kernel
    primeTerm_theta(d) = sum_{n = p^k} (log p) n^{-1/2} (n^{1/4} + n^{-1/4})
                                       (triG(log n + d) + triG(-log n + d)).
EXACT (elementary): n^{1/4} + n^{-1/4} >= 2 for all n >= 1 (AM-GM) and
triG >= 0 on R, both kernels being sums of nonnegative terms, hence
primeTerm_theta(d) >= 2 * primeTerm(d) for every lag d. The control defect is
    D_theta(c) := archQ(c) - kform(primeTerm_theta)(c)
                  - 2 * bigPsiI2^2 * (sum_i c_i sinh(x_i/2))^2,
i.e. ONLY the prime kernel is substituted; the archimedean kernel and the
pole correction stay at zeta's. D_theta is a numerical control object: it is
NOT the full Weil form of F_theta (whose own archimedean kernel and pole
structure are the separate named gap above) and NOT a Lean Prop.

KILL CRITERION 2 (adopted 2026-11-12): the apparatus is validated iff the
F_theta control goes NEGATIVE at the windows where zeta's margins are
smallest (the clustered L = 0.5 family; the near-coalescent Nelder-Mead
minimiser; the n = 2 family). If D_theta stays >= 0 everywhere zeta survives,
the machinery's silence on zeta means nothing and the ArchDom numerics are
indicted, not SemilocalArchDomAll.

Runs:
  0. cross-checks: fast vs slow theta kernel; domination check; n=2 weight.
  1. n=2 family a in {0.05,0.3,1,2,4,5,6,8,10}, c = (cosh(a/2), -1), L = a:
     D_zeta and D_theta side by side.
  2. projected-eigenvalue scans: scout grid L in {0.5,1,2,4} and extension
     grid L in {5,8} (configs equispaced / clustered0.1L / endpoints, n=3..6):
     lam_min of D_theta vs D_zeta on ker l_c.
  3. the extension's near-coalescent Nelder-Mead point (n=4, read from its
     results JSON): lam_min of D_theta vs D_zeta.
  4. mpmath dps=40 confirmation: n=2 at a in {0.3,1,4}; projected D_theta
     spectrum for clustered0.1L L=0.5 n=6 (zeta's smallest confirmed case).
Verdict threshold (scout convention): negative means < -5e-6 (float),
with dps40 confirmation of the load-bearing cases.
"""
import json
import math
import sys
import time

import numpy as np

sys.path.insert(0, ".")
import research.flint_archdom_scout_2026_11_12 as scout  # noqa: E402
import research.flint_archdom_margin_extension_2026_11_12 as ext  # noqa: E402

BIGPSI_I2 = scout.BIGPSI_I2
TH = 0.25
OUT = {"label": "NUMERICAL", "experiment": "flint_ftheta_control",
       "datum": "F_theta = xi(s-1/4)xi(s+1/4); prime weights Lambda(n) n^{-1/2} (n^{1/4}+n^{-1/4}); "
                "archimedean kernel and pole correction kept at zeta's (record-pinned scope)"}


# ---------------- theta prime kernel ----------------
def prime_kernel_theta_fast(d):
    d = float(d)
    lim = math.exp(2 + abs(d))
    s = 0.0
    for n, lp in ext.PP:
        if n > lim:
            break
        ln = math.log(n)
        w = n ** TH + n ** (-TH)
        s += lp / math.sqrt(n) * w * (scout.triG(ln + d) + scout.triG(-ln + d))
    return s


def prime_kernel_theta_slow(d):
    """Direct factoring version, scout.prime_kernel with the theta weights."""
    d = float(d)
    nmax = int(math.exp(2 + abs(d))) + 2
    s = 0.0
    for n in range(2, nmax + 1):
        p, _ = scout._factor_prime(n)
        if p is None:
            continue
        ln = math.log(n)
        w = n ** TH + n ** (-TH)
        s += math.log(p) / math.sqrt(n) * w * (scout.triG(ln + d) + scout.triG(-ln + d))
    return s


_theta_cache = {}
def theta_lag(d):
    k = round(abs(float(d)), 10)
    if k not in _theta_cache:
        _theta_cache[k] = prime_kernel_theta_fast(k)
    return _theta_cache[k]


# ---------------- wide prime-power table (lags up to 14) ----------------
def prime_power_table16(nmax):
    logp = {}
    sieve = bytearray(b"\x01") * (nmax + 1)
    sieve[0] = 0
    sieve[1] = 0
    for p in range(2, nmax + 1):
        if not sieve[p]:
            continue
        lp = math.log(p)
        q = p
        while q <= nmax:
            logp[q] = lp
            if q > nmax // p:
                break
            q *= p
        if p * p <= nmax:
            start = p * p
            sieve[start:nmax + 1:p] = b"\x00" * ((nmax - start) // p + 1)
    return sorted(logp.items())


PP16 = prime_power_table16(int(math.exp(16)) + 3)  # covers lags |d| <= 14


def prime_kernel_fast16(d):
    d = float(d)
    lim = math.exp(2 + abs(d))
    s = 0.0
    for n, lp in PP16:
        if n > lim:
            break
        ln = math.log(n)
        s += lp / math.sqrt(n) * (scout.triG(ln + d) + scout.triG(-ln + d))
    return s


def prime_kernel_theta_fast16(d):
    d = float(d)
    lim = math.exp(2 + abs(d))
    s = 0.0
    for n, lp in PP16:
        if n > lim:
            break
        ln = math.log(n)
        w = n ** TH + n ** (-TH)
        s += lp / math.sqrt(n) * w * (scout.triG(ln + d) + scout.triG(-ln + d))
    return s


def build_D_theta(x):
    x = np.asarray(x, float)
    n = len(x)
    A = np.zeros((n, n)); Pth = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            d = abs(x[i] - x[j])
            A[i, j] = scout.lag_val(scout._arch_cache, lambda k: scout.arch_scipy(k)[0], round(d, 10))
            Pth[i, j] = theta_lag(d)
    s = np.sinh(x / 2)
    Dth = A - Pth - 2 * BIGPSI_I2 ** 2 * np.outer(s, s)
    return A, Dth


def crosschecks():
    print("== Part 0: theta kernel cross-checks ==")
    w2 = 2 ** TH + 2 ** (-TH)
    print(f"  weight at n=2: {w2:.16f} (banked 2.0301035302564356...)")
    OUT["weight_n2"] = w2
    worst = 0.0
    for d in [0.0, 0.3, 1.0, 2.0, 4.0, 6.0]:
        a, b = prime_kernel_theta_slow(d), prime_kernel_theta_fast(d)
        z = ext.prime_kernel_fast(d)
        worst = max(worst, abs(a - b))
        print(f"  d={d}: theta slow={a:+.12f} fast={b:+.12f} |diff|={abs(a-b):.3e}; "
              f"theta/zeta ratio={a/z:.6f} (domination >= 2: {a >= 2*z})")
        assert a >= 2 * z - 1e-12, "domination primeTerm_theta >= 2 primeTerm violated"
    OUT["theta_fast_vs_slow_maxdiff"] = worst
    print(f"  fast vs slow max |diff| = {worst:.3e}")
    assert worst < 1e-9
    # wide-table (e^16) agreement: zeta vs scout slow (small d) and vs e^14 table (d <= 12)
    for d in [0.0, 0.3, 1.0, 2.0, 4.0]:
        assert abs(prime_kernel_fast16(d) - scout.prime_kernel(d)) < 1e-9
        assert abs(prime_kernel_theta_fast16(d) - prime_kernel_theta_slow(d)) < 1e-9
    for d in [6.0, 10.0, 12.0]:
        assert abs(prime_kernel_fast16(d) - ext.prime_kernel_fast(d)) < 1e-9
        assert abs(prime_kernel_theta_fast16(d) - prime_kernel_theta_fast(d)) < 1e-9
    print("  wide-table e^16 kernels agree with scout slow + e^14 table where ranges overlap: OK")


def part1():
    print("\n== Part 1: n=2 family, D_zeta vs D_theta at c = (cosh(a/2), -1) ==")
    print(f"{'a':>5} {'archQ':>13} {'primeQ_z':>13} {'primeQ_th':>13} {'D_zeta':>11} {'D_theta':>11}")
    rows = []
    for a in [0.05, 0.3, 1.0, 2.0, 4.0, 5.0, 6.0, 8.0, 10.0]:
        x = np.array([0.0, a])
        c = np.array([np.cosh(a / 2), -1.0])
        A, PmS, Pm, Dz, K, md = scout.build_mats(x, a)
        _, Dth = build_D_theta(x)
        archQ = float(c @ A @ c)
        pq_z = float(c @ PmS @ c)
        s = np.sinh(x / 2)
        polecorr = 2 * BIGPSI_I2 ** 2 * float(c @ s) ** 2
        pq_th = archQ - polecorr - float(c @ Dth @ c)
        Dv_z = float(c @ Dz @ c)
        Dv_th = float(c @ Dth @ c)
        # FS1 consistency (oracle advice 2026-11-13, section A): D - D_theta must
        # equal -(primeQ_theta - primeQ_zeta) exactly (arch + pole shared).
        fs1 = abs((Dv_z - Dv_th) - (pq_th - pq_z))
        print(f"{a:>5} {archQ:>+13.8f} {pq_z:>+13.8f} {pq_th:>+13.8f} {Dv_z:>+11.8f} {Dv_th:>+11.8f}  (FS1 resid {fs1:.2e})")
        assert fs1 < 1e-9, "FS1 identity violated: control build has arch/pole leakage"
        rows.append({"a": a, "archQ": archQ, "primeQ_zeta": pq_z, "primeQ_theta": pq_th,
                     "pole_corr": polecorr, "D_zeta": Dv_z, "D_theta": Dv_th, "fs1_resid": fs1,
                     "verdict_theta": "NEGATIVE (control detects)" if Dv_th < -5e-6 else "non-negative"})
    OUT["n2"] = rows


def part1b():
    """FS2 (oracle advice 2026-11-13, section B; ORACLE SIGHTING, unaudited):
    the n=2 family at a = log p, with the lag-0 / lag-log p decomposition.
    Sighting's prediction: D_zeta > 0 throughout; D_theta < 0 for large p.
    Either half failing kills the mechanism. mpmath dps=40 is primary for
    p >= 1009; scipy values are reported alongside."""
    print("\n== Part 1b: FS2 family x = (0, log p), lag decomposition ==")
    print(f"{'p':>7} {'a=log p':>8} {'D_zeta':>11} {'D_theta':>11} {'lag0_z':>11} {'laga_z':>11} "
          f"{'lag0_th':>11} {'laga_th':>11} {'polecorr':>10} {'spike_z':>10} {'spike_th':>10}")
    small_primes = [p for p in range(2, 98) if scout._factor_prime(p)[0] == p]
    rows = []
    for p in small_primes + [1009, 10007, 100003, 999983]:
        a = math.log(p)
        x = np.array([0.0, a])
        c = np.array([np.cosh(a / 2), -1.0])
        if p <= 97:
            A, PmS, Pm, Dz, K, md = scout.build_mats(x, a)
            _, Dth = build_D_theta(x)
            Pth_m = np.array([[theta_lag(abs(x[i] - x[j])) for j in range(2)] for i in range(2)])
        else:
            # scout's factoring kernel is too slow at lag log p for large p;
            # build directly from the wide-table kernels (cross-checked in part 0)
            A00, Aa = scout.arch_scipy(0.0)[0], scout.arch_scipy(a)[0]
            A = np.array([[A00, Aa], [Aa, A00]])
            pz0, pza = prime_kernel_fast16(0.0), prime_kernel_fast16(a)
            pt0, pta = prime_kernel_theta_fast16(0.0), prime_kernel_theta_fast16(a)
            PmS = np.array([[pz0, pza], [pza, pz0]])
            Pth_m = np.array([[pt0, pta], [pta, pt0]])
            poleRank = 2 * BIGPSI_I2 ** 2 * np.outer(np.sinh(x / 2), np.sinh(x / 2))
            Dz = A - PmS - poleRank
            Dth = A - Pth_m - poleRank
        s = np.sinh(x / 2)
        polecorr = 2 * BIGPSI_I2 ** 2 * float(c @ s) ** 2
        cc = float(c[0] ** 2 + c[1] ** 2)
        Mz = A - PmS
        Mth = A - Pth_m
        lag0_z = cc * Mz[0, 0]; laga_z = 2 * c[0] * c[1] * Mz[0, 1]
        lag0_t = cc * Mth[0, 0]; laga_t = 2 * c[0] * c[1] * Mth[0, 1]
        Dv_z = lag0_z + laga_z - polecorr
        Dv_th = lag0_t + laga_t - polecorr
        assert abs(Dv_z - float(c @ Dz @ c)) < 1e-9 and abs(Dv_th - float(c @ Dth @ c)) < 1e-9
        # the n = p spike inside primeTerm(a) / primeTerm_theta(a): triG(2 log p) + triG(0)
        tg = scout.triG(2 * a) + scout.triG(0.0)
        spike_z = math.log(p) / math.sqrt(p) * tg
        spike_t = math.log(p) / math.sqrt(p) * (p ** TH + p ** (-TH)) * tg
        print(f"{p:>7} {a:>8.4f} {Dv_z:>+11.8f} {Dv_th:>+11.8f} {lag0_z:>+11.8f} {laga_z:>+11.8f} "
              f"{lag0_t:>+11.8f} {laga_t:>+11.8f} {polecorr:>10.6f} {spike_z:>10.6f} {spike_t:>10.6f}")
        rows.append({"p": p, "a": a, "D_zeta": Dv_z, "D_theta": Dv_th,
                     "lag0_zeta": lag0_z, "laga_zeta": laga_z,
                     "lag0_theta": lag0_t, "laga_theta": laga_t,
                     "pole_corr": polecorr, "spike_zeta_at_p": spike_z,
                     "spike_theta_at_p": spike_t})
    OUT["fs2_family_scipy"] = rows


def scan_grid(Ls, tag):
    print(f"\n== {tag}: projected lam_min on ker l_c, D_zeta vs D_theta ==")
    print(f"{'L':>4} {'config':>14} {'n':>2} {'lam_zeta':>13} {'lam_theta':>13} {'mu_theta':>10}")
    rows = []
    for L in Ls:
        for n in range(3, 7):
            for name, x in scout.configs(L, n).items():
                A, PmS, Pm, Dz, K, md = scout.build_mats(x, L)
                _, Dth = build_D_theta(x)
                lam_z, _, _ = scout.proj_eigmin(Dz, x)
                lam_t, _, _ = scout.proj_eigmin(Dth, x)
                normA = float(np.linalg.norm(A, 2))
                print(f"{L:>4} {name:>14} {n:>2} {lam_z:>+13.9f} {lam_t:>+13.9f} {lam_t/normA:>+10.3e}")
                rows.append({"L": L, "config": name, "n": n, "lam_min_zeta": float(lam_z),
                             "lam_min_theta": float(lam_t), "mu_theta": float(lam_t / normA),
                             "verdict_theta": "NEGATIVE (control detects)" if lam_t < -5e-6 else "non-negative"})
    OUT[tag] = rows


def part3():
    print("\n== Part 3: near-coalescent Nelder-Mead point from the margin extension ==")
    res = json.load(open("research/flint_archdom_margin_extension_2026_11_12_results.json"))
    row = min(res["minimisation_n4"], key=lambda r: r["lam_min"])
    x = np.array(row["x"], float)
    L = float(row["L"])
    A, PmS, Pm, Dz, K, md = scout.build_mats(x, L)
    _, Dth = build_D_theta(x)
    lam_z, _, _ = scout.proj_eigmin(Dz, x)
    lam_t, _, _ = scout.proj_eigmin(Dth, x)
    print(f"  x={np.round(x, 4)} L={L}: lam_zeta={lam_z:+.9f} lam_theta={lam_t:+.9f}")
    OUT["nm_point"] = {"x": [float(v) for v in x], "L": L,
                       "lam_min_zeta": float(lam_z), "lam_min_theta": float(lam_t),
                       "verdict_theta": "NEGATIVE (control detects)" if lam_t < -5e-6 else "non-negative"}


def part4():
    print("\n== Part 4: mpmath confirmation (dps=40) ==")
    import mpmath as mp
    mp.mp.dps = 40

    def triG_mp(y):
        ay = abs(y)
        if ay <= 1:
            return mp.mpf(2) / 3 - y * y + ay ** 3 / 2
        if ay <= 2:
            return (2 - ay) ** 3 / 6
        return mp.mpf(0)

    PPMP = [(n, mp.log(int(round(math.exp(lp))))) for n, lp in PP16]

    def prime_mp(d, theta):
        d = mp.mpf(d)
        lim = mp.e ** (2 + abs(d))
        s = mp.mpf(0)
        for n, lp_mp in PPMP:
            if n > lim:
                break
            ln = mp.log(n)
            w = mp.mpf(1)
            if theta:
                w = mp.mpf(n) ** mp.mpf("0.25") + mp.mpf(n) ** mp.mpf("-0.25")
            s += lp_mp / mp.sqrt(n) * w * (triG_mp(ln + d) + triG_mp(-ln + d))
        return s

    def arch_mp(d):
        d = mp.mpf(d)
        f = lambda t: (mp.mpf(1) if t == 0 else (mp.sin(t / 2) / (t / 2)) ** 2) ** 2 \
            * mp.cos(t * d) * (mp.re(mp.digamma(mp.mpc(mp.mpf("0.25"), t / 2))) - mp.log(mp.pi))
        if d >= 2:
            nper = int(100 * d / mp.pi) + 1
            pts = [mp.mpf(0)] + [mp.mpf(k) * mp.pi / d for k in range(1, nper + 1)] + [mp.inf]
            return (1 / mp.pi) * mp.quad(f, pts)
        return (1 / mp.pi) * mp.quad(f, [0, mp.inf])

    bI2 = 16 * mp.sinh(mp.mpf("0.25")) ** 2
    conf = {"n2_mp": [], "cluster_mp": None}
    A0 = arch_mp("0")
    for a in ["0.3", "1", "4"]:
        av = mp.mpf(a)
        Aa = arch_mp(a)
        c0, c1 = mp.cosh(av / 2), mp.mpf(-1)
        archQ = c0 * c0 * A0 + 2 * c0 * c1 * Aa + A0
        polecorr = 2 * bI2 ** 2 * (c1 * mp.sinh(av / 2)) ** 2
        Dz = archQ - (prime_mp("0", False) * (c0 * c0 + 1) + 2 * c0 * c1 * prime_mp(a, False) + polecorr)
        Dt = archQ - (prime_mp("0", True) * (c0 * c0 + 1) + 2 * c0 * c1 * prime_mp(a, True) + polecorr)
        print(f"  a={a}: D_zeta={mp.nstr(Dz, 15)}  D_theta={mp.nstr(Dt, 15)}")
        conf["n2_mp"].append({"a": a, "D_zeta_mp": mp.nstr(Dz, 20), "D_theta_mp": mp.nstr(Dt, 20)})

    # projected D_theta spectrum, clustered0.1L L=0.5 n=6
    xlist = list(scout.configs(0.5, 6)["clustered0.1L"])
    xs = [mp.mpf(v) for v in xlist]
    n = len(xs)
    ell = [mp.cosh(v / 2) for v in xs]
    sv = [mp.sinh(v / 2) for v in xs]
    _amp, _pt = {}, {}
    def Amp(d):
        k = mp.nstr(d, 15)
        if k not in _amp:
            _amp[k] = arch_mp(d)
        return _amp[k]
    def Pth(d):
        k = mp.nstr(d, 15)
        if k not in _pt:
            _pt[k] = prime_mp(d, True)
        return _pt[k]
    Dm = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            d = abs(xs[i] - xs[j])
            Dm[i, j] = Amp(d) - Pth(d) - 2 * bI2 ** 2 * sv[i] * sv[j]
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
    Mp_ = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            Mp_[i, j] = (Z[i].T * Dm * Z[j])[0]
    ev = mp.eigsy(Mp_, eigvals_only=True)
    print(f"  clustered0.1L L=0.5 n=6: mp D_theta projected lam_min = {mp.nstr(min(ev), 15)} "
          f"(zeta float/mp anchor: +8.88e-7)")
    conf["cluster_mp"] = {"case": "clustered0.1L L=0.5 n=6", "lam_min_theta_mp": mp.nstr(min(ev), 20)}

    # FS2 large primes at dps40 (primary values for the sign prediction)
    fs2 = []
    for p in [1009, 10007, 100003, 999983]:
        av = mp.log(p)
        Aa = arch_mp(av)
        c0, c1 = mp.cosh(av / 2), mp.mpf(-1)
        archQ = c0 * c0 * A0 + 2 * c0 * c1 * Aa + A0
        polecorr = 2 * bI2 ** 2 * (c1 * mp.sinh(av / 2)) ** 2
        Dz = archQ - (prime_mp("0", False) * (c0 * c0 + 1) + 2 * c0 * c1 * prime_mp(av, False) + polecorr)
        Dt = archQ - (prime_mp("0", True) * (c0 * c0 + 1) + 2 * c0 * c1 * prime_mp(av, True) + polecorr)
        print(f"  FS2 p={p}: D_zeta={mp.nstr(Dz, 15)}  D_theta={mp.nstr(Dt, 15)}")
        fs2.append({"p": p, "D_zeta_mp": mp.nstr(Dz, 20), "D_theta_mp": mp.nstr(Dt, 20)})
    conf["fs2_mp"] = fs2
    OUT["mp_confirmation"] = conf


if __name__ == "__main__":
    t0 = time.time()
    crosschecks()
    part1()
    part1b()
    scan_grid([0.5, 1.0, 2.0, 4.0], "scan_scout_grid")
    scan_grid([5.0, 8.0], "scan_extension_grid")
    part3()
    part4()
    OUT["runtime_s"] = time.time() - t0
    # FS2 verdict (sighting prediction: D_zeta > 0 throughout, D_theta < 0 for large p;
    # mp values primary where present)
    fs2_mp = {r["p"]: r for r in OUT["mp_confirmation"]["fs2_mp"]}
    fs2_zeta_ok = all(r["D_zeta"] > 0 for r in OUT["fs2_family_scipy"]) and \
        all(float(r["D_zeta_mp"]) > 0 for r in fs2_mp.values())
    fs2_theta_neg_large = all(float(fs2_mp[p]["D_theta_mp"]) < 0 for p in [10007, 100003, 999983])
    OUT["fs2_verdict"] = {
        "zeta_half_holds": bool(fs2_zeta_ok),
        "theta_negative_at_large_p_mp": bool(fs2_theta_neg_large),
        "verdict": ("FS2 PREDICTION CONFIRMED (NUMERICAL): D_zeta > 0 on the whole log p family; "
                    "D_theta < 0 at p = 10007, 100003, 999983 (dps40) — the sighting's race "
                    "mechanism survives its named falsifier"
                    if fs2_zeta_ok and fs2_theta_neg_large else
                    "FS2 PREDICTION FAILS (NUMERICAL) in at least one half — see fs2 rows; "
                    "the race mechanism as stated is killed or needs its range corrected")}
    print(f"\nFS2: {OUT['fs2_verdict']['verdict']}")
    neg_windows = []
    for key in ["scan_scout_grid", "scan_extension_grid"]:
        neg_windows += [f"{r['config']} L={r['L']} n={r['n']}" for r in OUT.get(key, [])
                        if r["verdict_theta"].startswith("NEGATIVE")]
    neg_n2 = [r["a"] for r in OUT["n2"] if r["verdict_theta"].startswith("NEGATIVE")]
    detects_at_smallest = any(r["config"] == "clustered0.1L" and r["L"] == 0.5
                              and r["verdict_theta"].startswith("NEGATIVE")
                              for r in OUT.get("scan_scout_grid", []))
    OUT["kill_criterion_2"] = {
        "detects_at_clustered_L0.5": bool(detects_at_smallest),
        "negative_scan_windows": neg_windows,
        "negative_n2_a": neg_n2,
        "verdict": ("CONTROL PASSES (NUMERICAL): the apparatus detects the F_theta defect "
                    "at zeta's smallest-margin windows"
                    if detects_at_smallest else
                    "CONTROL FAILS (NUMERICAL): D_theta does not go negative where zeta's "
                    "margins are smallest; the apparatus is not discriminating at those windows")}
    print(f"\nKILL CRITERION 2: {OUT['kill_criterion_2']['verdict']}  (runtime {OUT['runtime_s']:.1f}s)")
    with open("research/flint_ftheta_control_2026_11_13_results.json", "w") as fh:
        json.dump(OUT, fh, indent=1)
