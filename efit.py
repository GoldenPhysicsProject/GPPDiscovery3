#!/usr/bin/env python3
"""E(c): how well does the CvS ground-state transform approximate Xi?

This is the quantity in Connes 2026 section 6 (open problem). CvS Thm 6.1
already gives real-rootedness of the finite characteristic function for every
c; what is missing is that these converge to Xi. Hurwitz then finishes.

    E(c) = sup_tau |F_even(tau)/F_even(0) - Xi(tau)/Xi(0)| / max|Xi/Xi(0)|

Measured, not proved. If E(c) decays with a clean shape, that shape is the
candidate theorem. Needs no spectral gap, no eigenvector tracking, no prolate
candidate, no Hilbert space: it is a comparison of two functions each computed
to 60+ digits.

At c=7, N=28 we measured E = 4.84e-2 with only the primes 2,3,5,7, and the
error did NOT grow as the tau-window widened from 15 to 30 to 60. That
flatness is the signature of locally uniform convergence.
"""
import argparse, json, time
import mpmath as mp
import connes_cvs as cc

p = argparse.ArgumentParser()
p.add_argument("--c", type=int, required=True)
p.add_argument("--N", type=int, required=True)
p.add_argument("--T", type=int, default=300)
p.add_argument("--dps", type=int, default=60)
p.add_argument("--taumax", type=float, default=60.0)
p.add_argument("--out", default="result.json")
a = p.parse_args()

t0 = time.time()
mp.mp.dps = a.dps
Q = cc.build_galerkin_matrix(a.c, N=a.N, T=a.T, dps=a.dps)
lam, v = cc.compute_ground_state(Q)
L = mp.log(a.c)
N = a.N

def F_even(tau):
    tot = mp.mpc(0)
    for k in range(-N, N + 1):
        den = 2 * mp.pi * k / L - tau
        gk = L if abs(den) < mp.mpf('1e-30') else (mp.e**(-1j * tau * L) - 1) / (1j * den)
        tot += v[k + N, 0] * gk
    return mp.re(mp.e**(1j * tau * L / 2) * tot) / mp.sqrt(L)

def xi(s):
    s = mp.mpc(s)
    return s * (s - 1) * mp.pi**(-s / 2) * mp.gamma(s / 2) * mp.zeta(s)

taus = [mp.mpf(k) / 4 for k in range(0, int(4 * a.taumax) + 1)]
F = [F_even(t) for t in taus]
X = [mp.re(xi(mp.mpf('0.5') + 1j * t)) for t in taus]
Fn = [f / F[0] for f in F]
Xn = [x / X[0] for x in X]
scale = max(abs(x) for x in Xn)

def E_upto(tmax):
    idx = [i for i, t in enumerate(taus) if t <= tmax]
    return float(max(abs(Fn[i] - Xn[i]) for i in idx) / scale)

rec = dict(c=a.c, N=a.N, T=a.T, dps=a.dps, a=float(L),
           lam1=mp.nstr(lam, 20), lam1_log10=float(mp.log(abs(lam), 10)),
           E_15=E_upto(15), E_30=E_upto(30), E_60=E_upto(a.taumax),
           F_at_g1=float(Fn[[i for i,t in enumerate(taus) if abs(t-mp.mpf('14.25'))<mp.mpf('0.13')][0]]),
           seconds=round(time.time() - t0, 1))
json.dump(rec, open(a.out, "w"))
print(json.dumps(rec), flush=True)
