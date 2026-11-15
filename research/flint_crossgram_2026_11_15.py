#!/usr/bin/env python
"""Flint eleventh run (2026-11-15): the two-centre cross Gram.

C2-a (run 9) consumes off-diagonal concentration Grams
W[ (1/2) A_e(rho,.) A_e(rho',.) ] for rho != rho', whose spectral
face is cross(a, b, e) := sum_w m_w Ahat^a(w) Ahat^b(w) with
a = tau rho, b = tau rho'.  Heart's T2 is the same tail estimate
for the centre pair (q, tau q).  This script computes
cross(rho_1, rho_2, e) at e = 0.12 and 0.06 two ways:
  (i) direct zero sum over staged Platt ordinates |gamma| <= 3000;
  (ii) residues: R = Ahat^a Ahat^b has simple poles at a+-e, b+-e
       with Res(Ahat^q, q+e) = -e/2, Res(Ahat^q, q-e) = +e/2
       (Heart J2 principal parts), so by J1
       cross = -sum_p Res(R L, p)
             = (e/2)[ Ahat^b(a+e) L(a+e) - Ahat^b(a-e) L(a-e)
                    + Ahat^a(b+e) L(b+e) - Ahat^a(b-e) L(b-e) ].
The third leg -- run 9's arithmetic assembly C_e(rho_1, rho_2)/2 --
is computed in flint_mevaluator_compare_2026_11_15.py (addendum).
Prediction (agreed TailLemma, T2 two-centre form): cross -> 0 as
e -> 0 at rate O(e^2 log) at fixed separation |rho_1 - rho_2|.

Epistemic status: NUMERICAL (dps 50).
"""

import json
import sys

sys.path.insert(0, '/home/hatch/workspace/zeros-platt')

import mpmath
from mpmath import mp, mpf, mpc, log, digamma

mp.dps = 50
from platt_sequential import zeros_in_range, zeros_by_index

GAM = sorted(zeros_in_range(0, 3000))
G12 = dict(zeros_by_index(1, 2))
RHO1 = mpc(mpf('0.5'), mpf(G12[1]))
RHO2 = mpc(mpf('0.5'), mpf(G12[2]))
ZXI = []
for g in GAM:
    ZXI.append(mpc(mpf('0.5'), g))
    ZXI.append(mpc(mpf('0.5'), -g))


def Ahat(q, w, e):
    return e**2 / (e**2 - (w - q)**2)


def L_xi(w):
    return (mp.zeta(w, derivative=1) / mp.zeta(w) + 1 / w + 1 / (w - 1)
            - log(mp.pi) / 2 + digamma(w / 2) / 2)


def cross_residue(a, b, e):
    return (e / 2) * (Ahat(b, a + e, e) * L_xi(a + e)
                      - Ahat(b, a - e, e) * L_xi(a - e)
                      + Ahat(a, b + e, e) * L_xi(b + e)
                      - Ahat(a, b - e, e) * L_xi(b - e))


OUT = {'rho1': mpmath.nstr(RHO1, 20), 'rho2': mpmath.nstr(RHO2, 20), 'rows': []}
for e in (mpf('0.12'), mpf('0.06')):
    s = mpc(0)
    for w in ZXI:
        s += Ahat(RHO1, w, e) * Ahat(RHO2, w, e)
    r = cross_residue(RHO1, RHO2, e)
    row = {'eps': str(e), 'cross_zerosum': mpmath.nstr(s, 20),
           'cross_residue': mpmath.nstr(r, 20),
           'diff': mpmath.nstr(abs(s - r), 4)}
    OUT['rows'].append(row)
    print('cross(rho1, rho2) eps', e, ' sum', row['cross_zerosum'],
          ' residue', row['cross_residue'], ' diff', row['diff'])

with open('/home/hatch/workspace/gppdiscovery3/research/'
          'flint_crossgram_2026_11_15_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('done')
