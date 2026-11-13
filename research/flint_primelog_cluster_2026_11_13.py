#!/usr/bin/env python3
"""
NUMERICAL run of the prime-log cluster family against the open Prop
`SemilocalArchDomAll` (Flint / Discovery 3, 2026-11-13; queue item 4(i),
named as the next falsifier by flint_archdom_margin_extension).

Labels: NUMERICAL only. Statements verbatim in flint_archdom_scout_2026_11_12.py
(GPPVerify WeilGramContraction.lean / WeilGramArchPD.lean @ 5fb6c8d0):
  SemilocalArchDomAll := forall L > 0, SemilocalArchDom (primesBelowExp L) L,
  i.e. on ker l_c (sum_i c_i cosh(x_i/2) = 0), pairwise |x_i - x_j| <= L,
  D(c) := archQ x c - primePoleQ (primesBelowExp L) x c >= 0.
On ker l_c, D(c) = kform Kpsi x c exactly (PROVED, kform_KS_of_ker).
A single negative projected eigenvalue beyond error KILLS SemilocalArchDomAll
(and, by the banked hPsi_iff_semilocalArchDomAll, H_psi) at that instance
(Kill Criterion 1). No finite positive evidence proves the Prop.

WHY THIS FAMILY. The defect form is NOT translation invariant: archTerm and
primeTerm see only the lags |x_i - x_j|, but the constraint ker l_c
(sum c_i cosh(x_i/2) = 0) and the pole correction (sum c_i sinh(x_i/2))^2 see
the absolute position. A cluster centred at u0 = log P therefore probes a
genuinely different slice of the statement than the origin-centred scans:
at large u0 the constraint weights cosh(u0/2) ~ sqrt(P)/2 and the pole
correction is suppressed by e^{-u0} on ker l_c. Prime-log points are also
where the prime kernel's own mass sits (triG(log n +/- d) samples), so a
cluster whose internal lags are prime-log differences is the configuration
class most aligned with the prime side.

FAMILIES (L is always the actual spread, S = primesBelowExp L, matched):
  A tight cluster at log P: x_i = log P + delta_i, delta equispaced in
    [-eps, eps], eps in {0.05, 0.125, 0.25}, n in {4, 6},
    P in {101, 1009, 10007, 100003, 1000003} (primes).
  B consecutive primes: x_i = log p_i for the first n primes >= P
    (n in {4, 6, 8}), P in {1009, 10007, 100003}; L = log p_n - log p_1.
  C prime-power offsets: x = log P + {-log 3, -log 2, 0, log 2, log 3},
    P in {1009, 100003}; L = 2 log 3. Internal lags are exactly log(p/q).

For every configuration the worst c on ker l_c is computed EXACTLY as the
projected min eigenvalue of the D-form (scout.proj_eigmin); the margin of
record mu(x) = lam_min / ||A(x)||_2 (spectral norm) is reported alongside.
Verdict threshold (scout convention): KILLED if lam_min < -5e-6 (float);
anything in [-1e-3, 0] is re-run at mpmath dps=40 before any verdict.
mpmath dps=40 confirmation is run for the five smallest-lam_min cases.
"""
import json
import math
import sys
import time

import numpy as np

sys.path.insert(0, ".")
import research.flint_archdom_scout_2026_11_12 as scout  # noqa: E402
import research.flint_archdom_margin_extension_2026_11_12 as ext  # noqa: E402

OUT = {"label": "NUMERICAL", "experiment": "flint_primelog_cluster", "rows": []}


def primes_upto(N):
    sieve = [True] * (N + 1)
    out = []
    for i in range(2, N + 1):
        if sieve[i]:
            out.append(i)
            if i * i <= N:
                for j in range(i * i, N + 1, i):
                    sieve[j] = False
    return out


def first_primes_from(P, n):
    """First n primes >= P (simple deterministic search)."""
    def is_prime(m):
        if m < 2:
            return False
        d = 2
        while d * d <= m:
            if m % d == 0:
                return False
            d += 1 if d == 2 else 2
        return True
    out, m = [], P
    while len(out) < n:
        if is_prime(m):
            out.append(m)
        m += 1
    return out


def run_config(name, x, L):
    x = np.asarray(x, float)
    A, PmS, Pm, D, K, md = scout.build_mats(x, L)
    lam, lamabs, cvec = scout.proj_eigmin(D, x)
    lamK, _, _ = scout.proj_eigmin(K, x)
    normA = float(np.linalg.norm(A, 2))
    row = {"case": name, "n": int(len(x)), "L": float(L),
           "x": [float(v) for v in x],
           "lam_min_D": float(lam), "lam_max_abs_D": float(lamabs),
           "lam_min_Kform": float(lamK), "norm_A": normA,
           "mu_margin_of_record": float(lam / normA) if normA > 0 else None,
           "primeS_vs_prime_maxdiff": float(md),
           "verdict": "SURVIVES (NUMERICAL)" if lam >= -5e-6 else "KILLED at this instance (NUMERICAL)"}
    OUT["rows"].append(row)
    print(f"{name:>44} n={len(x)} L={L:.4f} lam_min={lam:+.9f} |D|max={lamabs:.6f} "
          f"mu={row['mu_margin_of_record']:+.3e} (K-form {lamK:+.9f})")
    return row


def families():
    print("== Family A: tight clusters at log P ==")
    for P in [101, 1009, 10007, 100003, 1000003]:
        u0 = math.log(P)
        for eps in [0.05, 0.125, 0.25]:
            for n in [4, 6]:
                d = np.linspace(-eps, eps, n)
                run_config(f"A P={P} eps={eps} n={n}", u0 + d, 2 * eps)
    print("\n== Family B: consecutive-prime log clusters ==")
    for P in [1009, 10007, 100003]:
        for n in [4, 6, 8]:
            ps = first_primes_from(P, n)
            x = np.log(np.array(ps, float))
            run_config(f"B P={P} n={n} primes={ps[0]}..{ps[-1]}", x, float(x[-1] - x[0]))
    print("\n== Family C: prime-power log offsets ==")
    offs = np.array([-math.log(3), -math.log(2), 0.0, math.log(2), math.log(3)])
    for P in [1009, 100003]:
        run_config(f"C P={P} offsets=log(1/3..3)", math.log(P) + offs, 2 * math.log(3))


def mp_confirm(rows, topk=5):
    """mpmath dps=40 projected spectra for the smallest-lam_min cases."""
    print(f"\n== mpmath confirmation (dps=40): {topk} smallest lam_min cases ==")
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
        for n, lp in ext.PP:
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

    _amp, _pmp = {}, {}
    def Amp(d):
        k = mp.nstr(mp.mpf(d), 15)
        if k not in _amp:
            _amp[k] = arch_mp(mp.mpf(d))
        return _amp[k]
    def Pmp(d):
        k = mp.nstr(mp.mpf(d), 15)
        if k not in _pmp:
            _pmp[k] = prime_mp(mp.mpf(d))
        return _pmp[k]

    bI2 = 16 * mp.sinh(mp.mpf("0.25")) ** 2
    conf = []
    for row in sorted(rows, key=lambda r: r["lam_min_D"])[:topk]:
        xs = [mp.mpf(v) for v in row["x"]]
        n = len(xs)
        ell = [mp.cosh(v / 2) for v in xs]
        sv = [mp.sinh(v / 2) for v in xs]
        Dm = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                d = abs(xs[i] - xs[j])
                Dm[i, j] = Amp(d) - Pmp(d) - 2 * bI2 ** 2 * sv[i] * sv[j]
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
        mn = min(ev)
        print(f"{row['case']}: mp lam_min = {mp.nstr(mn, 15)}  (float {row['lam_min_D']:+.9f})")
        conf.append({"case": row["case"], "lam_min_mp": mp.nstr(mn, 20),
                     "lam_min_float": row["lam_min_D"],
                     "verdict_mp": "SURVIVES (NUMERICAL)" if mn >= 0 else "NEGATIVE at dps40 (NUMERICAL)"})
    OUT["mp_confirmation"] = conf


if __name__ == "__main__":
    t0 = time.time()
    families()
    mp_confirm(OUT["rows"])
    OUT["runtime_s"] = time.time() - t0
    OUT["primeS_vs_prime_maxdiff_overall"] = max(r["primeS_vs_prime_maxdiff"] for r in OUT["rows"])
    killed = [r for r in OUT["rows"] if r["verdict"].startswith("KILLED")]
    OUT["overall"] = ("KILLED (NUMERICAL) at a tested instance: " + "; ".join(r["case"] for r in killed)
                      if killed else "SURVIVES every tested instance (NUMERICAL)")
    print(f"\nOVERALL: {OUT['overall']}  (runtime {OUT['runtime_s']:.1f}s)")
    with open("research/flint_primelog_cluster_2026_11_13_results.json", "w") as fh:
        json.dump(OUT, fh, indent=1)
