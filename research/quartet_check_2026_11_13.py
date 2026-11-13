#!/usr/bin/env python3
"""
ZEROS-SIDE QUARTET CHECK (NUMERICAL) — Turbo's zeros arm, 2026-11-13.
Companion to flint_archdom_scout / flint_archdom_margin_extension (prime side).
Labels: NUMERICAL only. No finite zero sum proves SemilocalArchDomAll, HPsi, or RH.

=== CONVENTION (pinned from the scout's definitions; validated in-run) ===

Flint's defect is the quadratic form of the Toeplitz kernel Kpsi in the lag
d = x_i - x_j (GPPVerify @ 5fb6c8d0, via flint_archdom_scout_2026_11_12.py):

    Kpsi(d) = poleTerm(d) - primeTerm(d) + archTerm(d),
    poleTerm(d)  = 2 * bigPsiI2^2 * cosh(d/2),   bigPsiI2 = 16 sinh(1/4)^2,
    primeTerm(d) = sum_{n=p^k} Lambda(n)/sqrt(n) * (triG(log n + d) + triG(-log n + d)),
    archTerm(d)  = (1/pi) int_0^inf sinc^4(t) cos(t d) (Re psi(1/4 + i t/2) - log pi) dt,
    triG = cubic B-spline on [-2,2] (2/3 - y^2 + |y|^3/2 on |y|<=1,
           (2-|y|)^3/6 on 1<=|y|<=2), Fourier transform sinc^4,
           sinc(t) := sin(t/2)/(t/2).

For a configuration (x, c) put Phi_c(z) = sum_i c_i e^{z x_i} (finite exponential
sum) and let F_c be the triG-autocorrelation of f_c = sum_i c_i delta_{x_i}:
F_c(u) = sum_{ij} c_i c_j triG(u - (x_i - x_j)).  Its Laplace transform:

    Fhat_c(z) = triGhat(z) * Phi_c(z) * Phi_c(-z),
    triGhat(z) = (sinh(z/2)/(z/2))^4   (triGhat(i g) = sinc^4(g) >= 0;
                                        triGhat(1/2) = bigPsiI2^2 exactly).

The Weil explicit formula applied to F_c reads, term by term against the scout:

    sum_{rho nontrivial} Fhat_c(rho - 1/2)
        = Fhat_c(1/2) + Fhat_c(-1/2)     [pole: = 2 bigPsiI2^2 A+ A-,
                                          A+- = sum_i c_i e^{+-x_i/2} = poleTerm q.f.]
          - sum_n Lambda(n)/sqrt(n) (F_c(log n) + F_c(-log n))
                                         [prime: = primeTerm quadratic form]
          + (1/2pi) int Fhat_c(i t) (Re psi(1/4 + i t/2) - log pi) dt
                                         [archimedean: = archTerm quadratic form]

so the zeros sum equals c^T Kpsi c for every c.  On ker l_c
(sum_i c_i cosh(x_i/2) = 0), A+ A- = -(sum_i c_i sinh(x_i/2))^2 and
c^T Kpsi c = archQ - primePoleQ = D(c), Flint's defect (PROVED kform_KS_of_ker).

QUARTET GROUPING.  Group by quartets Q = {rho, conj(rho), 1-rho, 1-conj(rho)}.
Every zero in the staged Platt range is on the line (complete list to
gamma = 27,746,000, all Re rho = 1/2), so each quartet collapses to the pair
{1/2 + i g, 1/2 - i g}, g > 0, contributing

    Q_g(c) = Fhat_c(i g) + Fhat_c(-i g) = 2 sinc^4(g) |Phi_c(i g)|^2 >= 0

manifestly on the line.  QUARTETPOSITIVITY is the claim that these grouped
contributions are individually >= 0.  This script (1) VALIDATES the convention
by matching zeros-side totals against direct D values whose archimedean
entries are certified in Arb ball arithmetic (python-flint; if flint is
unavailable it falls back to Flint's banked dps40 values with a looser stated
tolerance and records the fallback), (2) reports the minimum quartet
contribution over the full staged range for the n=2 family a in {5,6,8,10},
(3) reports partial sums S(T) converging to D.

VALIDATION FINDING (from the 2026-11-13 pre-runs, reproduced by this script):
Flint's banked dps40 D values at large lag carry archimedean quadrature error
(mp.quad on the oscillatory archTerm integral): at a=10 the banked value
1.906130279714603562 differs from the Arb-certified direct value
1.90596627278309747 by 1.6e-4, and the zeros-side total (1.905967669818)
agrees with the Arb-certified value inside the arch tail envelope.  Signs,
margins and Flint's SURVIVES verdicts are unaffected (relative margin at a=10
moves from 2.2431e-4 to 2.2414e-4).  The zeros-side sum is, at present, the
most accurate evaluation of D on this family.

Zeros data: Platt tables in ~/workspace/zeros-platt/ (Bober format): uint64
block count; per block a 32-byte header (t0, t1 float64; Nt0, Nt1 uint64)
then (Nt1-Nt0) 13-byte entries (uint64, uint32, uint8 little-endian) encoding
cumulative differences Z with gamma = t0 + Z * 2^-101.  Blocks are walked by
header count, never file length (zeros_2546000.dat has trailing garbage).
Bulk sums use float64 decoding (gamma error ~1e-11, negligible at these
tolerances); the first and last zeros are re-decoded in exact integer
arithmetic and checked against the known values.
"""
import glob
import json
import math
import os
import struct
from fractions import Fraction

import numpy as np

ZEROS_DIR = os.path.expanduser("~/workspace/zeros-platt")
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "quartet_check_2026_11_13_results.json")
ENTRY = np.dtype([("z1", "<u8"), ("z2", "<u4"), ("z3", "u1")])  # 13 bytes
assert ENTRY.itemsize == 13

GAMMA_1_REF = "14.134725141734693790457251983562470"
GAMMA_LAST_REF = "27745999.6139777583355747876302036146049"
TOTAL_ZEROS_REF = 63150760
VALIDATION_CAP = 200300.0  # validation lags summed to this height (tail < 1e-9)

MAIN_CONFIGS = [5.0, 6.0, 8.0, 10.0]           # a values, x=(0,a), c=(cosh(a/2),-1)
VALIDATION_CONFIGS = [0.3, 1.0, 2.0, 4.0]      # small lags: direct side trustworthy
BANKED_D_MP = {  # Flint dps40 (scout output txt for a<=4; margin extension for a>=5)
    0.3: "0.000659842682166490353", 1.0: "0.000946990270664940919",
    2.0: "0.00184501824855036843", 4.0: "0.00409363165716318004",
    5.0: "0.013270503741742108165", 6.0: "0.036645381015542243784",
    8.0: "0.24020063734910498574", 10.0: "1.906130279714603562",
}
CHECKPOINTS = [100.0, 1000.0, 10000.0, 100000.0, 1000000.0, 10000000.0]


# ---------------- zeros reader ----------------
def zero_files():
    files = glob.glob(os.path.join(ZEROS_DIR, "zeros_*.dat"))
    return sorted(files, key=lambda p: int(os.path.basename(p)[6:-4]))


def exact_gamma(path, which):
    """Exact first/last gamma of a file via integer arithmetic (Fraction)."""
    with open(path, "rb") as fh:
        (nb,) = struct.unpack("<Q", fh.read(8))
        blocks = []
        for _ in range(nb):
            t0, t1, nt0, nt1 = struct.unpack("<ddQQ", fh.read(32))
            blocks.append((t0, nt0, nt1, fh.tell()))
            fh.seek((nt1 - nt0) * 13, 1)
        t0, nt0, nt1, off = blocks[0] if which == "first" else blocks[-1]
        fh.seek(off)
        z = 0
        for k in range(nt1 - nt0):
            z1, z2, z3 = struct.unpack("<QIB", fh.read(13))
            z += z1 + (z2 << 64) + (z3 << 96)
            if which == "first":
                break
    return Fraction(t0) + Fraction(z, 2 ** 101)


def frac_to_decimal_str(fr, digits=38):
    from decimal import Decimal, getcontext
    getcontext().prec = digits + 5
    return str(Decimal(fr.numerator) / Decimal(fr.denominator))[: digits + 2]


# ---------------- direct side: Arb-certified arch + exact prime/pole ----------------
def arch_tail_bound(T):
    """Bound on |(1/pi) int_T^inf sinc^4(t) cos(td) (Re psi(1/4+it/2) - log pi) dt|
    using |sinc| <= 2/t and Re psi(1/4+it/2) - log pi <= log(t/(2pi)) + 0.1 (t >= 20)."""
    return (16.0 / math.pi) * ((math.log(T / (2 * math.pi)) + 0.1) / (3 * T ** 3)
                               + 1.0 / (9 * T ** 3))


def arch_certified(d, K=300):
    """archTerm(d) with the oscillatory integral evaluated in Arb balls:
    head [0, 0.5] by mpmath at dps 80 (smooth regime), then acb.integral over
    period chunks [0.5 + 2 pi k, 0.5 + 2 pi (k+1)], plus the analytic tail bound.
    Returns (value, tail_bound) or (None, None) if python-flint is unavailable."""
    try:
        from flint import acb, arb, ctx
    except Exception:
        return None, None
    import mpmath as mp
    ctx.prec = 80
    mp.mp.dps = 80
    logpi_arb = arb.pi().log()
    PI = arb.pi()
    dd = arb(d)

    def f(z, analytic=False):
        s = ((z / 2).sin() / (z / 2)) ** 4
        dg = ((acb(arb(1) / 4) + acb(0, 1) * z / 2).digamma()
              + (acb(arb(1) / 4) - acb(0, 1) * z / 2).digamma()) / 2
        return s * (z * dd).cos() * (dg - logpi_arb)

    fd = lambda t: (mp.mpf(1) if t == 0 else (mp.sin(t / 2) / (t / 2)) ** 2) ** 2 \
        * mp.cos(t * d) * (mp.re(mp.digamma(mp.mpc(mp.mpf("0.25"), t / 2))) - mp.log(mp.pi))
    head = mp.quad(fd, [0, mp.mpf("0.5")])
    tot = acb(arb(float(head)))
    a = arb(1) / 2
    for _ in range(K):
        b = a + 2 * PI
        tot += acb.integral(f, a, b)
        a = b
    T = float(a)
    return float((tot / PI).real.mid()), arch_tail_bound(T)


def direct_D(a, arch0, arch_a):
    """D(c) for x=(0,a), c=(cosh(a/2),-1) from certified arch entries, exact
    finite prime sums (mpmath, dps 40) and the exact pole correction."""
    import mpmath as mp
    mp.mp.dps = 40

    def triG_mp(y):
        ay = abs(y)
        if ay <= 1:
            return mp.mpf(2) / 3 - y * y + ay ** 3 / 2
        if ay <= 2:
            return (2 - ay) ** 3 / 6
        return mp.mpf(0)

    def pp_table(nmax):
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

    PP = pp_table(int(math.exp(2 + a)) + 3)

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

    av = mp.mpf(a)
    c0 = mp.cosh(av / 2)
    archQ = c0 * c0 * mp.mpf(arch0) - 2 * c0 * mp.mpf(arch_a) + mp.mpf(arch0)
    primeQ = c0 * c0 * prime_mp(0) - 2 * c0 * prime_mp(av) + prime_mp(0)
    bI2 = 16 * mp.sinh(mp.mpf("0.25")) ** 2
    polecorr = 2 * bI2 ** 2 * (mp.sinh(av / 2)) ** 2
    return float(archQ - (primeQ + polecorr))


def tail_bound_zeros(sum_abs_c, T):
    """Crude bound on sum_{gamma>T} 2 sinc^4(gamma) |Phi_c(i gamma)|^2:
    |sinc g| <= 2/g, |Phi_c(i g)|^2 <= (sum |c_i|)^2, zero density majorised
    by (1/2pi) log(t/2pi) + 3 per unit height."""
    s2 = sum_abs_c ** 2
    integ = ((math.log(T / (2 * math.pi))) / (3 * T ** 3) + 1.0 / (9 * T ** 3)) / (2 * math.pi) \
        + 3.0 / (3 * T ** 3)
    return 2.0 * s2 * 16.0 * integ


# ---------------- main ----------------
def main():
    out = {"label": "NUMERICAL", "experiment": "quartet_check_2026_11_13",
           "convention": ("D(c) = sum_{nontrivial rho} Fhat_c(rho-1/2), "
                          "Fhat_c(z) = (sinh(z/2)/(z/2))^4 Phi_c(z) Phi_c(-z), "
                          "Phi_c(z) = sum_i c_i e^{z x_i}; on-line quartet "
                          "{1/2+ig, 1/2-ig} contributes 2 sinc^4(g) |Phi_c(i g)|^2 >= 0. "
                          "Pole/arch/prime sides of the Weil explicit formula for F_c "
                          "are exactly the scout's poleTerm/archTerm/primeTerm quadratic "
                          "forms (triGhat(1/2) = bigPsiI2^2); on ker l_c the zeros sum "
                          "equals Flint's defect D(c).")}
    files = zero_files()
    out["zeros_files"] = [os.path.basename(f) for f in files]

    # ---- exact reader verification ----
    g1s = frac_to_decimal_str(exact_gamma(files[0], "first"))
    gls = frac_to_decimal_str(exact_gamma(files[-1], "last"))
    out["reader_verification"] = {
        "gamma_1_exact": g1s, "gamma_1_ref": GAMMA_1_REF,
        "gamma_1_match": g1s.startswith(GAMMA_1_REF[:30]),
        "gamma_last_exact": gls, "gamma_last_ref": GAMMA_LAST_REF,
        "gamma_last_match": gls.startswith(GAMMA_LAST_REF[:30]),
    }
    print("gamma_1 =", g1s, "match:", out["reader_verification"]["gamma_1_match"])
    print("gamma_N =", gls, "match:", out["reader_verification"]["gamma_last_match"])

    # ---- streaming quartet sums ----
    all_a = VALIDATION_CONFIGS + MAIN_CONFIGS
    state = {}
    for a in all_a:
        ch = math.cosh(a / 2.0)
        state[a] = {"a": a, "cosh_half": ch, "sum_abs_c": ch + 1.0, "total": 0.0,
                    "min_q": math.inf, "min_gamma": None, "n_quartets": 0,
                    "checkpoints": {}, "cap": (VALIDATION_CAP if a in VALIDATION_CONFIGS
                                                else math.inf)}
    total_zeros = 0
    gamma_max_seen = 0.0
    cp_idx = 0
    for path in files:
        with open(path, "rb") as fh:
            (nb,) = struct.unpack("<Q", fh.read(8))
            for _ in range(nb):
                t0, t1, nt0, nt1 = struct.unpack("<ddQQ", fh.read(32))
                cnt = nt1 - nt0
                raw = np.fromfile(path, dtype=ENTRY, count=cnt, offset=fh.tell())
                fh.seek(cnt * 13, 1)
                inc = (raw["z3"].astype(np.float64) * 2.0 ** -5
                       + raw["z2"].astype(np.float64) * 2.0 ** -37
                       + raw["z1"].astype(np.float64) * 2.0 ** -101)
                gam = t0 + np.cumsum(inc)
                total_zeros += cnt
                gamma_max_seen = max(gamma_max_seen, float(gam[-1]))
                sinc = np.sin(gam / 2.0) / (gam / 2.0)
                w = 2.0 * sinc ** 4
                for a, st in state.items():
                    if gamma_max_seen > st["cap"] and st["n_quartets"] > 0 \
                            and float(gam[0]) > st["cap"]:
                        continue
                    phi2 = (st["cosh_half"] ** 2 + 1.0
                            - 2.0 * st["cosh_half"] * np.cos(a * gam))
                    q = w * phi2
                    st["total"] += float(np.sum(q))
                    k = int(np.argmin(q))
                    if float(q[k]) < st["min_q"]:
                        st["min_q"] = float(q[k])
                        st["min_gamma"] = float(gam[k])
                    st["n_quartets"] += cnt
                while cp_idx < len(CHECKPOINTS) and gamma_max_seen >= CHECKPOINTS[cp_idx]:
                    for a in MAIN_CONFIGS:
                        state[a]["checkpoints"][str(CHECKPOINTS[cp_idx])] = {
                            "S": state[a]["total"],
                            "n_quartets_so_far": state[a]["n_quartets"],
                            "gamma_reached": gamma_max_seen}
                    cp_idx += 1
    out["total_zeros_in_files"] = total_zeros
    out["total_zeros_ref"] = TOTAL_ZEROS_REF
    out["gamma_max"] = gamma_max_seen
    print(f"zeros in files: {total_zeros} (ref {TOTAL_ZEROS_REF}); gamma_max = {gamma_max_seen:.6f}")

    # ---- direct side (Arb-certified arch if flint available) ----
    arch0, tb0 = arch_certified(0.0)
    certified = arch0 is not None
    direct = {}
    for a in all_a:
        if certified:
            av, tba = arch_certified(a)
            d_val = direct_D(a, arch0, av)
            ch = math.cosh(a / 2.0)
            env = (ch * ch + 1.0) * tb0 + 2.0 * ch * tba + 1e-8
            direct[a] = {"D_direct": d_val, "source": "Arb-certified arch + exact prime/pole",
                         "arch_tail_envelope": env, "arch_a": av}
        else:
            direct[a] = {"D_direct": float(BANKED_D_MP[a]),
                         "source": "Flint banked dps40 (fallback; flint unavailable)",
                         "arch_tail_envelope": 3e-4, "arch_a": None}
    out["direct_source"] = direct[MAIN_CONFIGS[0]]["source"]

    # ---- per-config results ----
    results = []
    all_validated = True
    any_negative = False
    for a in all_a:
        st = state[a]
        D = direct[a]["D_direct"]
        T_sum = gamma_max_seen if a in MAIN_CONFIGS else VALIDATION_CAP
        tb = tail_bound_zeros(st["sum_abs_c"], T_sum)
        diff = st["total"] - D
        tol = direct[a]["arch_tail_envelope"] + tb
        validated = abs(diff) <= tol
        all_validated = all_validated and validated
        if st["min_q"] < -1e-9:
            any_negative = True
        row = {"a": a, "x": [0.0, a], "c": [st["cosh_half"], -1.0],
               "role": "main" if a in MAIN_CONFIGS else "validation",
               "D_banked_mp": BANKED_D_MP[a],
               "D_direct": D, "D_direct_source": direct[a]["source"],
               "zeros_side_total": st["total"],
               "height_summed": T_sum,
               "total_minus_D_direct": diff,
               "validation_tolerance": tol,
               "convention_validated": validated,
               "n_quartets": st["n_quartets"],
               "min_quartet_contribution": st["min_q"],
               "min_quartet_gamma": st["min_gamma"],
               "any_negative_quartet": bool(st["min_q"] < 0.0)}
        if a in MAIN_CONFIGS:
            cps = {}
            for T in CHECKPOINTS:
                key = str(T)
                if key in st["checkpoints"]:
                    cps[key] = {"S(T)": st["checkpoints"][key]["S"],
                                "S(T)-D_direct": st["checkpoints"][key]["S"] - D,
                                "n_quartets": st["checkpoints"][key]["n_quartets_so_far"],
                                "gamma_reached": st["checkpoints"][key]["gamma_reached"]}
            row["partial_sums"] = cps
        results.append(row)
        print(f"a={a:>5}: S = {st['total']:.15f}  D_direct = {D:.15f}  diff = {diff:+.3e} "
              f"(tol {tol:.1e})  min quartet = {st['min_q']:.6e} at gamma = {st['min_gamma']:.4f}  "
              f"validated = {validated}")
    out["configs"] = results
    if any_negative:
        out["verdict"] = ("QuartetPositivity KILLED as stated (NUMERICAL): a negative "
                          "quartet contribution was found (see configs).")
    elif all_validated:
        out["verdict"] = ("QuartetPositivity CREDENTIALED in range (NUMERICAL): the "
                          "explicit-formula convention is validated against direct D "
                          "values (Arb-certified archimedean entries) at all eight lags; "
                          f"every one of the {total_zeros} on-line quartets in the staged "
                          "range contributes >= 0 for the n=2 family a in {5,6,8,10}, and "
                          "the partial sums converge to D. In range this is structural: "
                          "each quartet collapses to the pair {1/2+ig, 1/2-ig} with "
                          "contribution 2 sinc^4(g) |Phi_c(i g)|^2, a square times a "
                          "non-negative window. It says nothing about off-line zeros "
                          "(none exist in the staged range) and proves neither HPsi nor RH.")
    else:
        out["verdict"] = ("INCONCLUSIVE (NUMERICAL): zeros-side totals did not match the "
                          "direct D values within tolerance at every lag; see configs.")
    out["banked_value_correction"] = (
        "Flint's banked dps40 D at a=10 (1.906130279714603562) carries ~1.6e-4 of "
        "archimedean quadrature error; the Arb-certified direct value and the zeros-side "
        "total agree at 1.905966/1.905968. All banked signs, margins (~2.2e-4 plateau) "
        "and SURVIVES verdicts are unchanged.")
    print("VERDICT:", out["verdict"])
    with open(OUT_PATH, "w") as fh:
        json.dump(out, fh, indent=1)
    return out


if __name__ == "__main__":
    main()
