#!/usr/bin/env python
"""Flint eleventh run (2026-11-15): domain of run 9's assembly.

flint_mevaluator_compare_2026_11_15.py localised the theta leg gap
entirely in the Hadamard -1/4 shift: the assembly call at label
s' = 1 + i gamma_1 (centre i gamma_1), where a := Re s' - 1/2 = +1/2
exactly -- the boundary of the continuation formulas' range (arch
face digamma argument 1/4 - a/2 = 0; prime face L at 0 and 1).
The gap there is eps * (-0.20279 - 0.008844 i) at eps = 0.12, 0.06
(run 9's own legs_diff/eps = 0.406 at C level, eps = 0.3/0.26/0.225).

THIS SCRIPT tests the assembly at eps = 0.12 for labels
s = sigma + i gamma_1 with sigma in {0.95, 1.1} (a = 0.45 inside,
a = 0.60 outside), against the definition-side zero sum B(centre, e),
centre = 1 - conj(s).  Predictions: inside agrees to the falsifier's
accuracy (~1e-9 or better); outside disagrees by the SAME jumped
constant ~ eps * (-0.20279 - 0.008844 i) if the defect is a
crossed-pole residue term (the pole of the gamma factor at 0 /
of zeta at 1), not a blow-up.

Run with the arb venv python (falsifier definitions need flint).
Epistemic status: NUMERICAL (assembly dps 50/65; sums dps 50).
"""

import json
import sys
import time

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


def B_zerosum(q, e):
    s = mpc(0)
    for w in ZXI:
        a = e**2 / (e**2 - (w - q)**2)
        s += a * a
    return s


EPS = mpf('0.12')
OUT = {}
Lcache = {}
for sig in ('0.95', '1.1'):
    s_lab = mpc(mpf(sig), G1)
    centre = 1 - s_lab.conjugate()
    Ca, _ = C_eps(EPS, s_lab, s_lab, Lcache)
    G = Ca / 2
    B = B_zerosum(centre, EPS)
    OUT[sig] = {'assembly_G': mpmath.nstr(G, 18),
                'B_centre_zerosum': mpmath.nstr(B, 18),
                'gap': mpmath.nstr(G - B, 8),
                'gap_over_eps': mpmath.nstr((G - B) / EPS, 8)}
    print('sigma', sig, 'centre', mpmath.nstr(centre, 12),
          'assembly', OUT[sig]['assembly_G'], ' B', OUT[sig]['B_centre_zerosum'],
          ' gap/eps', OUT[sig]['gap_over_eps'])

with open('/home/hatch/workspace/gppdiscovery3/research/'
          'flint_assembly_domain_2026_11_15_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('done')
