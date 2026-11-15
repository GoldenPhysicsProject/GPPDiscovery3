#!/usr/bin/env python
"""Flint eleventh run (2026-11-15): isolate the theta-assembly gap.

flint_mevaluator_compare found: run 9's Hadamard assembly for the
F_theta diagonal at s_theta differs from the zero-sum B(tau s_theta)
by -0.203 eps - 0.0088 eps i (eps = 0.12, 0.06), while the zeta leg
at rho_1 agrees to 1.5e-12.  The theta assembly is
  C^theta = C_e(s_fake(u_th - 1/4)) + C_e(s_fake(u_th + 1/4)),
  u_th = 1/4 + i gamma_1,  s_fake(u') = 1 - conj(u'),
i.e. centres u' = i gamma_1 (label s' = 1 + i gamma_1) and
u' = 1/2 + i gamma_1 (label rho_1).  The spectral-face identity is
  B_Fth(tau s_theta, e) = B_xi(i gamma_1, e) + B_xi(rho_1, e)
(shifting the F_theta zero sum by -+1/4).  This script computes
B_xi(i gamma_1, e) by direct zero sum and by J2 residues at
e = 0.12, 0.06, so the failing leg is identified exactly:
assembly part for centre i gamma_1 (from the compare run's JSON)
vs B_xi(i gamma_1, e).

Epistemic status: NUMERICAL (dps 50).
"""

import json
import sys

sys.path.insert(0, '/home/hatch/workspace/zeros-platt')

import mpmath
from mpmath import mp, mpf, mpc, log, digamma, polygamma

mp.dps = 50
from platt_sequential import zeros_in_range, zeros_by_index

GAM = sorted(zeros_in_range(0, 3000))
G1 = mpf(dict(zeros_by_index(1, 1))[1])
ZXI = []
for g in GAM:
    ZXI.append(mpc(mpf('0.5'), g))
    ZXI.append(mpc(mpf('0.5'), -g))
Q = mpc(0, G1)  # centre i gamma_1 (sigma = 0, not a zero)


def Ahat(q, w, e):
    return e**2 / (e**2 - (w - q)**2)


def L_xi(w):
    return (mp.zeta(w, derivative=1) / mp.zeta(w) + 1 / w + 1 / (w - 1)
            - log(mp.pi) / 2 + digamma(w / 2) / 2)


def Lp_xi(w):
    z = mp.zeta(w)
    zp = mp.zeta(w, derivative=1)
    zpp = mp.zeta(w, derivative=2)
    return (zpp / z - (zp / z) ** 2 - 1 / w**2 - 1 / (w - 1) ** 2
            + polygamma(1, w / 2) / 4)


def B_res(q, e):
    return -((e**2 / 4) * (Lp_xi(q + e) + Lp_xi(q - e))
             + (e / 4) * (L_xi(q - e) - L_xi(q + e)))


OUT = {}
for e in (mpf('0.12'), mpf('0.06')):
    s = mpc(0)
    for w in ZXI:
        a = Ahat(Q, w, e)
        s += a * a
    r = B_res(Q, e)
    OUT[str(e)] = {'B_zerosum': mpmath.nstr(s, 20),
                   'B_residue': mpmath.nstr(r, 20),
                   'diff': mpmath.nstr(abs(s - r), 4)}
    print('B_xi(i gamma_1) eps', e, 'sum', OUT[str(e)]['B_zerosum'],
          ' residue', OUT[str(e)]['B_residue'], ' diff', OUT[str(e)]['diff'])

with open('/home/hatch/workspace/gppdiscovery3/research/'
          'flint_b_i_gamma1_2026_11_15_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('done')
