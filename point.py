#!/usr/bin/env python3
"""One (c, N, T, dps) point of the truncated-Weil scan.

Emits lambda_1, lambda_2 and the SPECTRAL GAP, not just the ground state.

WHY THE GAP MATTERS
The live proof strategy (Davis-Kahan / residual) needs
    r_c / delta_c = O(c^{-1/2}),
where r_c is the residual of an explicit candidate against the localized Weil
operator and delta_c = lambda_2 - lambda_1 is the spectral gap. That target is
POLYNOMIAL, not exponential, because Hurwitz only needs locally uniform
convergence on |Im z| < 1/2 (the critical strip), and Paley-Wiener on support
log c costs e^{(log c) Y} = c^Y with Y < 1/2. So the gap is a load-bearing
quantity and lambda_1 alone is not enough.

Parity: the CvS operator commutes with k -> -k. The candidate ground state is
even, so the relevant gap is the EVEN-sector gap. We also report the odd sector
so we can see whether an odd eigenvalue falls inside it.
"""
import argparse, json, platform, time
import mpmath as mp
import connes_cvs as cc


def sector_spectrum(Q, parity):
    """Full spectrum of Q restricted to the even or odd sector."""
    DIM = Q.rows
    N = (DIM - 1) // 2
    inv_sqrt2 = 1 / mp.sqrt(2)
    if parity == "even":                       # e_0, (e_k + e_{-k})/sqrt2
        V = mp.matrix(DIM, N + 1)
        V[N, 0] = mp.mpf(1)
        for k in range(1, N + 1):
            V[N + k, k] = inv_sqrt2
            V[N - k, k] = inv_sqrt2
    else:                                      # (e_k - e_{-k})/sqrt2
        V = mp.matrix(DIM, N)
        for k in range(1, N + 1):
            V[N + k, k - 1] = inv_sqrt2
            V[N - k, k - 1] = -inv_sqrt2
    eigs, _ = mp.eigsy(V.T * Q * V)
    return sorted([eigs[i] for i in range(len(eigs))])


p = argparse.ArgumentParser()
p.add_argument("--c", type=int, required=True)
p.add_argument("--N", type=int, required=True)
p.add_argument("--T", type=int, default=300)
p.add_argument("--dps", type=int, default=90)
p.add_argument("--out", default="result.json")
a = p.parse_args()

t0 = time.time()
mp.mp.dps = a.dps
Q = cc.build_galerkin_matrix(a.c, N=a.N, T=a.T, dps=a.dps)
ev_e = sector_spectrum(Q, "even")
ev_o = sector_spectrum(Q, "odd")
l1, l2 = ev_e[0], ev_e[1]
gap = l2 - l1
l10 = lambda x: float(mp.log(abs(x), 10)) if x != 0 else None

rec = dict(
    c=a.c, N=a.N, T=a.T, dps=a.dps, dim=2 * a.N + 1, a=float(mp.log(a.c)),
    lam1=mp.nstr(l1, 25), lam2=mp.nstr(l2, 25), gap=mp.nstr(gap, 25),
    odd1=mp.nstr(ev_o[0], 25),
    lam1_log10=l10(l1), lam2_log10=l10(l2), gap_log10=l10(gap),
    odd1_log10=l10(ev_o[0]),
    sign1=int(mp.sign(l1)),
    odd_inside_gap=bool(l1 < ev_o[0] < l2),
    seconds=round(time.time() - t0, 1),
    python=platform.python_version(), mpmath=mp.__version__,
    connes_cvs=getattr(cc, "__version__", "?"),
)
json.dump(rec, open(a.out, "w"))
print(json.dumps(rec), flush=True)
