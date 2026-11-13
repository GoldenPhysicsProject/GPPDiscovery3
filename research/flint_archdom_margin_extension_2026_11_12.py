#!/usr/bin/env python3
"""
NUMERICAL discriminating tests of the open Prop `SemilocalArchDomAll`
(Flint / Discovery 3, 2026-11-12, follow-up to flint_archdom_scout).

Labels: NUMERICAL only. Statements verbatim in flint_archdom_scout_2026_11_12.py
(GPPVerify WeilGramContraction.lean / WeilGramArchPD.lean @ 5fb6c8d0):
  SemilocalArchDomAll := forall L > 0, SemilocalArchDom (primesBelowExp L) L,
  i.e. on ker l_c (sum_i c_i cosh(x_i/2) = 0), with pairwise |x_i - x_j| <= L,
  D(c) := archQ x c - primePoleQ (primesBelowExp L) x c >= 0.
On ker l_c, D(c) = kform Kpsi x c exactly (PROVED, kform_KS_of_ker).
A single negative D beyond numerical error KILLS SemilocalArchDomAll (and, by
the banked hPsi_iff_semilocalArchDomAll, H_psi) at that instance. No finite
amount of positive evidence proves it.

Motivation (cheapest falsifier first): the scout covered windows a, L <= 4 and
found D > 0 with relative margin thinning 7.5e-2 -> 1.8e-4. An external oracle
sighting (Gemini, 2026-11-12, UNAUDITED) proposed, among other things:
  (i)   stress windows beyond the scouted range, a in [5, 10];
  (ii)  direct constrained minimisation over (x, c in ker l_c) instead of
        config scans;
  (iii) x clustered near log p for large p.
This experiment runs (i) and (ii); (iii) is filed as a named next falsifier
(the prime-log cluster family) in the bridge notes. The sighting's
"sinc-distributed c" idea is NOT used: for fixed x the worst c on ker l_c is
found exactly by the projected eigenvalue, which dominates any c-family scan.

Parts:
  0. fast prime-power table, cross-checked against the scout's prime_kernel.
  1. n=2 family x=(0,a), c=(cosh(a/2),-1), L=a, a in {5,6,8,10} (scipy route).
  2. projected min-eigenvalue scan at L in {5,8}, scout config families, n=3..6.
  3. Nelder-Mead minimisation over x of the projected min eigenvalue,
     n=4, L in {2,5}, two inits each (the direct falsifier).
  4. mpmath dps=40 confirmation of Part 1.
Verdict threshold: KILLED if any D or lam_min < -5e-6 (scout convention);
quadrature error estimates are printed with every scipy arch value used.
"""
import json
import math
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, ".")
import research.flint_archdom_scout_2026_11_12 as scout  # noqa: E402

BIGPSI_I2 = scout.BIGPSI_I2
OUT = {"label": "NUMERICAL", "experiment": "flint_archdom_margin_extension"}


# ---------------- fast prime kernel via a shared prime-power table ----------------
def prime_power_table(nmax):
    """[(n, log p)] for n = p^k <= nmax, by modified sieve (exact)."""
    logp = {}
    sieve = [True] * (nmax + 1)
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
            for j in range(p * p, nmax + 1, p):
                sieve[j] = False
    return sorted(logp.items())


PP = prime_power_table(int(math.exp(14)) + 3)  # covers lags |d| <= 12


def prime_kernel_fast(d):
    d = float(d)
    lim = math.exp(2 + abs(d))
    s = 0.0
    for n, lp in PP:
        if n > lim:
            break
        ln = math.log(n)
        s += lp / math.sqrt(n) * (scout.triG(ln + d) + scout.triG(-ln + d))
    return s


def crosschecks():
    print("== Part 0: fast prime kernel vs scout prime_kernel ==")
    worst = 0.0
    for d in [0.0, 0.3, 1.0, 2.0, 4.0, 6.0]:
        a, b = scout.prime_kernel(d), prime_kernel_fast(d)
        worst = max(worst, abs(a - b))
        print(f"  d={d}: scout={a:+.12f} fast={b:+.12f} |diff|={abs(a-b):.3e}")
    OUT["prime_fast_vs_scout_maxdiff"] = worst
    assert worst < 1e-9, "fast prime kernel disagrees with scout"
    a0 = scout.arch_scipy(0.0)
    print(f"  archTerm(0) scipy = {a0} (scout anchor: -1.542933... from the archpd files)")
    OUT["arch0_scipy"] = a0


# ---------------- Part 1: n=2 window extension ----------------
def n2_row(a, mats=None):
    x = np.array([0.0, a])
    c = np.array([np.cosh(a / 2), -1.0])
    A, PmS, Pm, D, K, md = mats if mats is not None else scout.build_mats(x, a)
    archQ = float(c @ A @ c)
    primeQ = float(c @ PmS @ c)
    polecorr = 2 * BIGPSI_I2 ** 2 * float(c @ np.sinh(x / 2)) ** 2
    P = primeQ + polecorr
    Dv = archQ - P
    Kv = float(c @ K @ c)
    return {"a": a, "archQ": archQ, "prime_form": primeQ, "pole_corr": polecorr,
            "primePoleQ": P, "D": Dv, "Kpsi_direct": Kv,
            "rel_margin_D_over_absA": Dv / abs(archQ),
            "rayleigh_D": Dv / float(c @ c),
            "primeS_vs_prime_maxdiff": md,
            "verdict": "SURVIVES (NUMERICAL)" if Dv >= -5e-6 else "KILLED at this instance (NUMERICAL)"}


def part1():
    print("\n== Part 1: n=2 family, window extension a in {5,6,8,10} (scipy) ==")
    print(f"{'a':>5} {'archQ':>13} {'P=primePoleQ':>13} {'D=A-P':>11} {'Kpsi_direct':>12} {'D/|A|':>10}")
    rows = []
    for a in [5.0, 6.0, 8.0, 10.0]:
        r = n2_row(a)
        rows.append(r)
        print(f"{a:>5} {r['archQ']:>+13.8f} {r['primePoleQ']:>+13.8f} {r['D']:>+11.8f} "
              f"{r['Kpsi_direct']:>+12.8f} {r['rel_margin_D_over_absA']:>+10.3e}")
    OUT["n2_extension_scipy"] = rows


# ---------------- Part 2: projected scan at L in {5,8} ----------------
def part2():
    print("\n== Part 2: projected min eigenvalue on ker l_c, L in {5,8} ==")
    print(f"{'L':>4} {'config':>14} {'n':>2} {'lam_min_D':>14} {'lam_max_abs_D':>14} {'rel_margin':>11}")
    rows = []
    for L in [5.0, 8.0]:
        for n in range(3, 7):
            for name, x in scout.configs(L, n).items():
                A, PmS, Pm, D, K, md = scout.build_mats(x, L)
                lam, lamabs, cvec = scout.proj_eigmin(D, x)
                print(f"{L:>4} {name:>14} {n:>2} {lam:>+14.9f} {lamabs:>14.6f} {lam/lamabs:>+11.3e}")
                rows.append({"L": L, "config": name, "n": n, "lam_min_D": float(lam),
                             "lam_max_abs_D": float(lamabs),
                             "verdict": "SURVIVES (NUMERICAL)" if lam >= -5e-6 else "KILLED at this instance (NUMERICAL)"})
    OUT["scan_L58"] = rows


# ---------------- Part 3: direct minimisation over x ----------------
def lam_min_of_x(x, L):
    x = np.asarray(x, float)
    spread = float(x.max() - x.min())
    A, PmS, Pm, D, K, md = scout.build_mats(x, L)
    lam, _, _ = scout.proj_eigmin(D, x)
    pen = 100.0 * max(0.0, spread - L) ** 2
    return lam + pen, lam, spread


def part3():
    print("\n== Part 3: Nelder-Mead minimisation of lam_min over x (n=4) ==")
    rows = []
    for L, inits in [(2.0, [np.linspace(-1, 1, 4), np.array([-0.1, -0.03, 0.03, 0.1])]),
                     (5.0, [np.linspace(-2.5, 2.5, 4), np.array([-0.25, -0.08, 0.08, 0.25])])]:
        for k, x0 in enumerate(inits):
            res = minimize(lambda v: lam_min_of_x(v, L)[0], x0, method="Nelder-Mead",
                           options={"maxiter": 160, "xatol": 1e-3, "fatol": 1e-7})
            _, lam, spread = lam_min_of_x(res.x, L)
            print(f"  L={L} init{k}: lam_min={lam:+.9f} at x={np.round(res.x, 4)} spread={spread:.4f} "
                  f"(nit={res.nit})")
            rows.append({"L": L, "init": k, "lam_min": float(lam),
                         "x": [float(v) for v in res.x], "spread": float(spread),
                         "verdict": "SURVIVES (NUMERICAL)" if lam >= -5e-6 else "KILLED at this instance (NUMERICAL)"})
    OUT["minimisation_n4"] = rows


# ---------------- Part 4: mpmath dps=40 confirmation of Part 1 ----------------
def part4():
    print("\n== Part 4: mpmath confirmation (dps=40), n=2 extension ==")
    import mpmath as mp
    mp.mp.dps = 40

    def triG_mp(y):
        ay = abs(y)
        if ay <= 1:
            return mp.mpf(2) / 3 - y * y + ay ** 3 / 2
        if ay <= 2:
            return (2 - ay) ** 3 / 6
        return mp.mpf(0)

    def prime_mp(d):
        d = mp.mpf(d)
        lim = mp.e ** (2 + abs(d))
        s = mp.mpf(0)
        for n, lp in PP:
            if n > lim:
                break
            ln = mp.log(n)
            s += mp.mpf(lp) / mp.sqrt(n) * (triG_mp(ln + d) + triG_mp(-ln + d))
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
    A0 = arch_mp("0")
    P0 = prime_mp("0")
    rows = []
    for a in ["5", "6", "8", "10"]:
        av = mp.mpf(a)
        Aa, Pa = arch_mp(a), prime_mp(a)
        c0, c1 = mp.cosh(av / 2), mp.mpf(-1)
        archQ = c0 * c0 * A0 + 2 * c0 * c1 * Aa + A0
        primeQ = c0 * c0 * P0 + 2 * c0 * c1 * Pa + P0
        polecorr = 2 * bI2 ** 2 * (c1 * mp.sinh(av / 2)) ** 2
        D = archQ - (primeQ + polecorr)
        rel = D / abs(archQ)
        print(f"  a={a}: D={mp.nstr(D, 18)}  D/|A|={mp.nstr(rel, 6)}  archQ={mp.nstr(archQ, 12)}")
        rows.append({"a": a, "D_mp": mp.nstr(D, 20), "rel_margin_mp": mp.nstr(rel, 8),
                     "archQ_mp": mp.nstr(archQ, 20)})
    OUT["n2_extension_mp"] = rows


if __name__ == "__main__":
    t0 = time.time()
    crosschecks()
    part1()
    part2()
    part3()
    part4()
    OUT["runtime_s"] = time.time() - t0
    verdicts = []
    for key in ["n2_extension_scipy", "scan_L58", "minimisation_n4"]:
        verdicts += [r["verdict"] for r in OUT.get(key, [])]
    OUT["overall"] = ("KILLED (NUMERICAL) at a tested instance"
                      if any(v.startswith("KILLED") for v in verdicts)
                      else "SURVIVES every tested instance (NUMERICAL)")
    print(f"\nOVERALL: {OUT['overall']}  (runtime {OUT['runtime_s']:.1f}s)")
    with open("research/flint_archdom_margin_extension_2026_11_12_results.json", "w") as fh:
        json.dump(OUT, fh, indent=1)
