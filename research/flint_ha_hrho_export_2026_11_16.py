#!/usr/bin/env python
"""Flint (2026-11-16): dyadic exports of run 16's (H-A) and (H-rho)
certificates, shaped to Claude's spec (GPP-bridge d60d24f) for his
Lean kernel checks, in the style of GPPVerify #606 / the item-4
dyadic export (bridge 00ef926).

(H-A): the 740 committed A_jj balls (j = 73..812) from run 16's
balls.json, hulled with the rebuilt-from-dumps balls, exported as
integer numerators over 2^160 (centre nearest, radius up). The
exported ball encloses the committed-dump ball (exact Fraction
check). c_A exported as the exact rational 867483/500000 = 1.734966
(<= run 16's certified 1.7349668350070715), plus run 16's committed
per-j profile lower bounds.

(H-rho): Atilde = D^{-1/2} A D^{-1/2} rebuilt in Arb balls from run
16's committed dumps via run 16's own verify path (shards + Q +
K2 balls -> assemble). Reports the max normalised off-diagonal
row-sum max_j sum_{k!=j} (|mid Atilde_jk| + rad Atilde_jk) and the
full Gershgorin bound including the diagonal deviation
|Atilde_jj - 1|. If that bound is < 1 it certifies rho_Z <= bound
(symmetric matrix: ||Atilde - I||_2 <= max Gershgorin deviation).
Exports the off-diagonal upper triangle on a 2^-32 grid (centres
nearest, radii up, all rounding folded into the radii) plus a
diag_dev array on the same grid, and re-derives the Gershgorin
bound from the exported integers alone (exact integer arithmetic).

Everything is computed from run 16's committed artifacts; no new
quadrature, no zero data. Epistemic status: the exports repackage
run 16's interval-certified balls; the Gershgorin step is an exact
consequence of those balls.
"""

import json
import math
import os
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
STAGE = '/tmp/flint_export'
os.makedirs(STAGE, exist_ok=True)

import importlib.util
spec = importlib.util.spec_from_file_location(
    'run16', os.path.join(HERE, 'flint_ha_hrho_cert_2026_11_15.py'))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

from flint import arb, arb_mat, ctx
ctx.prec = 200

T00 = time.time()
Z0, JMAX = 73, 812
N = JMAX - Z0 + 1
assert N == 740

# ------------------------------------------------------------ rebuild
B = m.load_shards()
Bph, Bpd, g, q_phi = m.shard_tables(B)
dump = json.load(open(m.BALLS_OUT))
Q = {int(k): arb(v) for k, v in dump['Q'].items()}
K2 = {int(k): arb(v) for k, v in json.load(open(m.K2_BALLS)).items()}
val = m.validate_shards(Bph, Bpd, g, Q)
assert val['V_all_ok'], 'run 16 validations failed on rebuild'
qd = m.q_delta_ball()
asm = m.assemble(Bph, Bpd, g, q_phi, qd, Q, K2)
A = asm['A']
print(f'rebuilt A from committed dumps ({time.time() - T00:.0f}s)',
      flush=True)

res = {'date': '2026-11-16', 'seat': 'Flint (D3)',
       'purpose': 'dyadic exports of run 16 (H-A)/(H-rho) for Claude '
                  'kernel checks, per bridge spec d60d24f',
       'validations_ok': True}

# ------------------------------------------------------------ (H-A)
GK_A = 160
GUARD = Fraction(1, 10 ** 36)          # covers ser() 30-digit parse slop
cen_num, rad_num, prof_str = [], [], []
max_diag_discrep = 0.0
for idx, j in enumerate(range(Z0, JMAX + 1)):
    ms, rs = dump['A_diag'][str(j)].split('+/-')
    mid_c, rad_c = Fraction(ms.strip()), Fraction(rs.strip())
    lo_c, hi_c = mid_c - rad_c, mid_c + rad_c
    rb = A[idx, idx]
    lo_r, hi_r = float(rb.lower()), float(rb.upper())
    max_diag_discrep = max(max_diag_discrep,
                           abs(float(rb.mid()) - float(mid_c)))
    lo = min(lo_c, Fraction(lo_r)) - GUARD
    hi = max(hi_c, Fraction(hi_r)) + GUARD
    cen = (lo + hi) / 2
    rad = (hi - lo) / 2
    cn = math.floor(cen * (1 << GK_A) + Fraction(1, 2))
    cerr = abs(cen - Fraction(cn, 1 << GK_A))
    rn = math.ceil((rad + cerr) * (1 << GK_A)) + 1
    # exact enclosure check against the committed ball
    assert Fraction(cn - rn, 1 << GK_A) <= lo_c
    assert Fraction(cn + rn, 1 << GK_A) >= hi_c
    cen_num.append(cn)
    rad_num.append(rn)
    prof_str.append(repr(dump['profile_lower'][str(j)]))
res['HA_diag_rebuild_max_mid_discrepancy'] = max_diag_discrep

# Arb profile check of the exported dyadic balls against c_A rational
CA_RAT = Fraction(867483, 500000)
min_prof = None
for idx, j in enumerate(range(Z0, JMAX + 1)):
    low = arb(cen_num[idx] - rad_num[idx]) / arb(1 << GK_A)
    w = m.om(j)
    p = low * w ** 4 / w.log()
    pl = float(p.lower())
    min_prof = pl if min_prof is None else min(min_prof, pl)
res['HA_export_profile_min'] = min_prof
assert min_prof >= float(CA_RAT), 'exported (H-A) balls fail c_A check'
print(f'(H-A) export done, profile min {min_prof} '
      f'({time.time() - T00:.0f}s)', flush=True)

ha = {'meta': {
        'quantity': 'A_jj balls, A = S_Z^eff of run 16 (L = 0.5), '
                    'j = 73..812',
        'source': 'GPPDiscovery3 flint_ha_hrho_cert_2026_11_15_balls.'
                  'json key A_diag (30-digit ser strings), hulled '
                  'with the same balls rebuilt from the committed '
                  'dumps; guard 1e-36 absolute',
        'format': 'value = num / 2^grid_k; ball = [cen-rad, cen+rad]; '
                  'centre rounded to nearest grid point, radius '
                  'rounded up so the exported ball encloses the '
                  'committed-dump ball (exact rational check)'},
      'grid_k': GK_A, 'j_range': [Z0, JMAX],
      'c_A_certified_run16': '1.7349668350070715',
      'c_A_rational': [867483, 500000],
      'A_cen_num': cen_num, 'A_rad_num': rad_num,
      'profile_lower_run16': prof_str}
json.dump(ha, open(os.path.join(STAGE, 'ha_dyadic_export_2026_11_16.json'),
                   'w'))
res['HA'] = {'grid_k': GK_A, 'entries': N,
             'c_A_rational': '867483/500000 = 1.734966',
             'profile_min_committed_run16': 1.7349668350070715,
             'profile_min_exact_from_committed_balls': min_prof,
             'discrepancy_note': 'run 16 committed profile_lower min '
             '(1.7349668350070715) sits 6.0e-9 BELOW the exact '
             'profile of run 16 committed A_diag balls '
             '(1.7349668410406225, verified at the argmin j=74 in '
             'mpmath at 80 digits). Conservative direction: the '
             'certified c_A understates the committed balls. No '
             'banked claim affected.'}

# ------------------------------------------------------------ (H-rho)
sq = [A[i, i].sqrt() for i in range(N)]
At = [[A[i, k] / (sq[i] * sq[k]) for k in range(N)] for i in range(N)]
print(f'Atilde built ({time.time() - T00:.0f}s)', flush=True)

# Symmetrise by pairwise hulls: the exact A is symmetric, but the
# Arb-solve balls differ by direction. Each hulled ball encloses
# both directional balls (hence the true entry); all row sums and
# the export use the hulled matrix.
Hs = [[None] * N for _ in range(N)]
max_asym = 0.0
for i in range(N):
    Hs[i][i] = At[i][i]
    for k in range(i + 1, N):
        max_asym = max(max_asym,
                       float(abs(At[i][k].mid() - At[k][i].mid()).upper())
                       + float(abs(At[i][k].rad() - At[k][i].rad()).upper()))
        h = At[i][k].union(At[k][i])
        Hs[i][k] = h
        Hs[k][i] = h
res['HRho_pair_ball_max_asymmetry'] = max_asym
At = Hs

# Arb row sums (off-diagonal) and diagonal deviations
rowsum = [arb(0)] * N
diagdev = [arb(0)] * N
for i in range(N):
    x = At[i][i]
    diagdev[i] = abs(x.mid() - 1) + x.rad()
    s = arb(0)
    for k in range(N):
        if k != i:
            x = At[i][k]
            s += abs(x.mid()) + x.rad()
    rowsum[i] = s
rs_up = [float(s.upper()) for s in rowsum]
dd_up = [float(d.upper()) for d in diagdev]
i_rs = max(range(N), key=lambda i: rs_up[i])
i_g = max(range(N), key=lambda i: rs_up[i] + dd_up[i])
res['HRho'] = {
    'offdiag_rowsum_max': rs_up[i_rs],
    'offdiag_rowsum_argmax_j': Z0 + i_rs,
    'diagdev_max': max(dd_up),
    'gershgorin_bound_arb': rs_up[i_g] + dd_up[i_g],
    'gershgorin_argmax_j': Z0 + i_g}
print(f'(H-rho) Arb row sums done: offdiag max {rs_up[i_rs]} at '
      f'j={Z0 + i_rs}, Gershgorin {rs_up[i_g] + dd_up[i_g]} '
      f'({time.time() - T00:.0f}s)', flush=True)

# Dyadic export on 2^-32 grid, upper triangle + diag_dev
GK_R = 32
GRID = 1 << GK_R
cen_l, rad_l = [], []
for i in range(N):
    for k in range(i + 1, N):
        x = At[i][k]
        cen_f = float(x.mid())
        delta = abs(x.mid() - arb(cen_f)) + x.rad()
        rad_frac = Fraction(float(delta.upper())) + Fraction(1, 2 ** 55)
        cen_frac = Fraction(cen_f)
        cn = math.floor(abs(cen_frac) * GRID + Fraction(1, 2))
        if cen_frac < 0:
            cn = -cn
        cerr = abs(cen_frac - Fraction(cn, GRID))
        rn = math.ceil((rad_frac + cerr) * GRID) + 1
        cen_l.append(cn)
        rad_l.append(rn)
dd_num = []
for i in range(N):
    dev = Fraction(dd_up[i]) + Fraction(1, 2 ** 55)
    dd_num.append(math.ceil(dev * GRID) + 1)
print(f'(H-rho) dyadic export built ({time.time() - T00:.0f}s)',
      flush=True)

# Exact integer self-check of the Gershgorin bound from the export
rows = [0] * N
pos = 0
for i in range(N):
    for k in range(i + 1, N):
        v = abs(cen_l[pos]) + rad_l[pos]
        rows[i] += v
        rows[k] += v
        pos += 1
G_num = max(dd_num[i] + rows[i] for i in range(N))
G_exp = G_num / GRID
res['HRho']['gershgorin_bound_from_export'] = G_exp
res['HRho']['gershgorin_bound_from_export_fraction'] = \
    f'{G_num}/{GRID}'
assert G_exp < 1.0, 'Gershgorin bound from export is vacuous'
print(f'(H-rho) export self-check: Gershgorin bound {G_exp} '
      f'({time.time() - T00:.0f}s)', flush=True)

hr = {'meta': {
        'quantity': 'Atilde = D^{-1/2} A D^{-1/2} of run 16, '
                    'D = diag(A); off-diagonal upper triangle, '
                    'indices 0..739 <-> j = 73..812, row-major '
                    '(i<k), plus diag_dev = upper bound of '
                    '|Atilde_jj - 1| per row',
        'source': 'A rebuilt in Arb balls from run 16 committed '
                  'dumps (R7 shards + Q balls + K2 balls) via run '
                  '16\'s own verify path (assemble)',
        'format': 'value = num / 2^32; ball = [cen-rad, cen+rad]; '
                  'centres are float64 midpoints of the Arb balls, '
                  'radii rounded up with every rounding (Arb->float, '
                  'float->grid) folded in',
        'certificate': 'rho_Z = ||Atilde - I||_2 <= max_j (diag_dev'
                       '_j + sum_{k!=j} (|cen_jk| + rad_jk)) = '
                       'gershgorin_bound_from_export (exact integer '
                       'arithmetic on this file)'},
      'grid_k': GK_R, 'dim': N, 'j_range': [Z0, JMAX],
      'offdiag_rowsum_max_arb': rs_up[i_rs],
      'gershgorin_bound_from_export': G_exp,
      'diag_dev_num': dd_num,
      'cen_num': cen_l, 'rad_num': rad_l}
json.dump(hr, open(os.path.join(
    STAGE, 'hrho_gershgorin_export_2026_11_16.json'), 'w'))
res['HRho']['grid_k'] = GK_R
res['HRho']['entries'] = len(cen_l)

# --------------------------------------- downstream room (NUMERICAL)
# Run 16 note section 4 components (filed floats):
NUM = 1.524894785662644e-11 + 4.907140201541864e-12
BETA0 = 1.30163e-8
RQ_D = (1 - 0.158) / 3.5942360966653345e-12   # from run 16 c_dir line
for tag, rho in (('at_export_gershgorin_bound', G_exp),
                 ('at_offdiag_rowsum_max', rs_up[i_rs])):
    v = NUM / (1 - rho)
    res.setdefault('downstream_NUMERICAL', {})[tag] = {
        'rho': rho, 'v0TLv0_le': v, 'room_vs_beta0': BETA0 / v,
        'c_dir_ge': (1 - rho) / RQ_D}
res['elapsed_s'] = time.time() - T00
json.dump(res, open(os.path.join(
    HERE, 'flint_ha_hrho_export_2026_11_16_results.json'), 'w'),
    indent=1)
print(json.dumps(res, indent=1))
print('EXPORT COMPLETE')
