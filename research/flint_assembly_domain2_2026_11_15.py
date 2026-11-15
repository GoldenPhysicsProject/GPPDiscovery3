#!/usr/bin/env python
"""Flint eleventh run (2026-11-15): pin the assembly's domain.

flint_assembly_domain found the assembly defect ~ -0.2028 eps at
centres i g1, 0.05 + i g1, -0.1 + i g1 (all with tau u != u), while
the on-line diagonal at rho_1 (tau u = u) agrees to 1.5e-12.
Hypothesis: the falsifier assembly's face formulas are derived for
centre pairs (u_k, u_j) with u_j = tau u_k (the pairing's shadow
symmetry); diagonal calls are then valid only at tau-fixed centres.

TESTS at eps = 0.12 (assembly vs definition-side zero sums):
 T1 diagonal, centre 1/2 + 20 i (tau-fixed, NOT a zero):
    assembly C_e(1/2+20i, 1/2+20i)/2  vs  B(1/2+20i, e).
 T2 tau-paired off-diagonal, centres (i g1, 1 + i g1):
    assembly C_e(1 + i g1, i g1)/2  vs  A_JGram(i g1, e)
    := sum_w m_w Ahat^{i g1}(w) Ahat^{1 + i g1}(w).
 T3 diagonal, centre 1/4 + i g1 = tau s_theta (not tau-fixed):
    assembly C_e(s_theta, s_theta)/2  vs  B(tau s_theta, e)
    = 1.0037359447046086 + 2.93727437461e-8 i (Script A).

Run with the arb venv python.  Epistemic status: NUMERICAL.
"""

import json
import sys

sys.path.insert(0, '/home/hatch/workspace/zeros-platt')

SRC = ('/home/hatch/workspace/gppdiscovery3/research/'
       'cesaro_gram_falsifier_2026_11_14.py')
text = open(SRC).read()
defs = text.split('# ================================================================ RUN')[0]
NS = {}
exec(compile(defs, SRC + ' [definitions]', 'exec'), NS)
mp = NS['mp']; mpf = NS['mpf']; mpc = NS['mpc']; mpmath = NS['mpmath']
C_eps = NS['C_eps']

from platt_sequential import zeros_in_range, zeros_by_index

GAM = sorted(zeros_in_range(0, 3000))
G1 = mpf(dict(zeros_by_index(1, 1))[1])
ZXI = []
for g in GAM:
    ZXI.append(mpc(mpf('0.5'), g))
    ZXI.append(mpc(mpf('0.5'), -g))


def Ahat(q, w, e):
    return e**2 / (e**2 - (w - q)**2)


def sumprod(a, b, e):
    s = mpc(0)
    for w in ZXI:
        s += Ahat(a, w, e) * Ahat(b, w, e)
    return s


EPS = mpf('0.12')
OUT = {}
Lcache = {}

# T1
lab = mpc(mpf('0.5'), mpf(20))
Ca, _ = C_eps(EPS, lab, lab, Lcache)
tgt = sumprod(lab, lab, EPS)
OUT['T1'] = {'assembly_G': mpmath.nstr(Ca / 2, 18),
             'target_B': mpmath.nstr(tgt, 18),
             'diff': mpmath.nstr(abs(Ca / 2 - tgt), 4)}
print('T1 tau-fixed non-zero centre:', OUT['T1'])

# T2: labels (1 + i g1, i g1) -> centres (i g1, 1 + i g1)
sk = mpc(mpf(1), G1)
sj = mpc(0, G1)
Ca, _ = C_eps(EPS, sk, sj, Lcache)
tgt = sumprod(mpc(0, G1), mpc(1, G1), EPS)
OUT['T2'] = {'assembly_G': mpmath.nstr(Ca / 2, 18),
             'target_A': mpmath.nstr(tgt, 18),
             'diff': mpmath.nstr(abs(Ca / 2 - tgt), 4)}
print('T2 tau-paired off-diagonal:', OUT['T2'])

# T3: label s_theta, centre tau s_theta = 1/4 + i g1
s_th = mpc(mpf('0.75'), G1)
Ca, _ = C_eps(EPS, s_th, s_th, Lcache)
tgt = mpc(mpf('1.003735944704608615663041'),
          mpf('2.937274374605970542614665e-8'))
OUT['T3'] = {'assembly_G': mpmath.nstr(Ca / 2, 18),
             'target_B_tau_stheta': mpmath.nstr(tgt, 18),
             'diff': mpmath.nstr(abs(Ca / 2 - tgt), 4)}
print('T3 diagonal at tau s_theta centre:', OUT['T3'])

with open('/home/hatch/workspace/gppdiscovery3/research/'
          'flint_assembly_domain2_2026_11_15_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('done')
