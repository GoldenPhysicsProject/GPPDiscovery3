#!/usr/bin/env python
"""Flint tenth run (2026-11-14): candidate (a) -- the shadow-antisymmetrised
concentration state in the COMPLETED Cesaro pairing.

QUESTION (Daniel's constructive turn; spec property P3 for candidate (a)).
Run 9 killed selection-by-extraction: the completed pairing's concentration
Gram at a zero tends to m_rho regardless of sigma (shadow-valued kernel:
the amplitude A_eps(s, w) is centred at w = 1 - conj(s)).  Candidate (a)
asks whether the combination ODD under the shadow exchange does better:
    v_odd = state(rho) - state(1 - rho),
the antisymmetric face of Daniel's boundary picture.  Hope: the quartet
spectral term becomes a square of a difference and the shadow mass cancels
by construction.  Test: compute the completed-pairing diagonal of v_odd
at zeta's first zero rho_1 (shadow partner 1 - rho_1 = conj(rho_1)) and at
F_theta's first zero s_th = rho_1 + 1/4 (partner 1 - s_th = 1/4 - i gamma_1),
plus the cross Gram G(rho, 1-rho) and the even combination for contrast.

PREDICTION (registered in rh_lab/catalog.json before the run): G_odd = C/2
tends to 2 at BOTH test zeros (the sum m_rho + m_{1-rho} of the two quartet
masses), cross Gram tends to 0, even combination has the same limit: NO
cancellation, NO selection -- the pairing's spectral measure is atomic on
the quartet, each atom fed by exactly one concentration centre, and a sign
between two states centred at different spectral points subtracts nothing
at either point.  If G_odd^theta <= 0 within bars while G_odd^zeta -> 2,
the prediction is wrong and candidate (a) stands.

EVALUATOR.  The falsifier arm's definitions, exec'd byte-for-byte (RUN
stripped), exactly as run 9 reused them.  C_eps(eps, s_k, s_j, Lcache)
assembles W[ (1/2) A(s_k,.) A(s_j,.) ] from pole - PRIME(continuation) +
ARCH.  F_theta pairings via Hadamard: W_theta[W_{ab}] assembled from the
two shifted centre pairs (u_a -/+ 1/4, u_b -/+ 1/4), as run 9 did for the
diagonal.  Spectral legs evaluate the combination weight
W_odd(w) = W_{aa}(w) - 2 W_{ab}(w) + W_{bb}(w) directly on the zero sets.

CERTIFICATION (before any odd number): C(0.3; rho_1, rho_1) arithmetic
against the falsifier's published 2.0000089211, and C_theta(a, a) at
eps = 0.30 against run 9's published rows (arith 2.51118750529894,
spectral 2.63282994493287).

Epistemic status: NUMERICAL throughout.
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
GAM = NS['GAM']; RHO = NS['RHO']
C_eps = NS['C_eps']; spectral_sum = NS['spectral_sum']
W_val = NS['W_val']; W_line = NS['W_line']; u_of = NS['u_of']
_quad = NS['quad']; _log = NS['log']

T_START = time.time()
OUT = {'evaluator': 'cesaro_gram_falsifier_2026_11_14.py definitions, exec\'d verbatim'}
Lcache = {}


def s_fake(u_prime):
    """Fake label whose amplitude centre u_of(s) equals u_prime (run 9)."""
    return 1 - u_prime.conjugate()


def C_theta(eps, a, b):
    """W_theta pair assembly: centres (u_a -/+ 1/4, u_b -/+ 1/4)."""
    ua, ub = u_of(a), u_of(b)
    tot = mpc(0)
    parts = {}
    for tag, da, db in (('-1/4', ua - mpf('0.25'), ub - mpf('0.25')),
                        ('+1/4', ua + mpf('0.25'), ub + mpf('0.25'))):
        v, _det = C_eps(eps, s_fake(da), s_fake(db), Lcache)
        parts[tag] = mpmath.nstr(v, 15)
        tot += v
    return tot, parts


def spectral_tail_bound(eps, T=3000):
    return _quad(lambda t: 2 * eps ** 4 * (_log(mpf(t)) + 13) / (mpf(t) - 30) ** 4,
                 [T, mpmath.inf])


def spec_odd_zeta(eps, ua, ub, T=3000):
    from platt_sequential import zeros_in_range
    S = mpc(0)
    for g in zeros_in_range(0, T):
        for sgn in (1, -1):
            t = sgn * g
            S += (W_line(t, eps, ua, ua) - 2 * W_line(t, eps, ua, ub)
                  + W_line(t, eps, ub, ub))
    return S


def spec_odd_theta(eps, ua, ub, T=3000):
    from platt_sequential import zeros_in_range
    S = mpc(0)
    for g in zeros_in_range(0, T):
        for re in (mpf('0.25'), mpf('0.75')):
            for sgn in (1, -1):
                w = mpc(re, sgn * g)
                S += (W_val(w, eps, ua, ua) - 2 * W_val(w, eps, ua, ub)
                      + W_val(w, eps, ub, ub))
    return S


# ------------------------------------------------------------ certification
print('== Certification of the reused evaluator ==')
C11_ar, _ = C_eps(mpf('0.3'), RHO[1], RHO[1], Lcache)
print('C(0.3; rho1, rho1) arithmetic =', mpmath.nstr(C11_ar, 18),
      ' published 2.0000089211')
OUT['cert_C11_arith'] = mpmath.nstr(C11_ar, 18)

a_th = RHO[1] + mpf('0.25')
Cth_aa_03, _ = C_theta(mpf('0.3'), a_th, a_th)
print('C_theta(0.3; a, a) arithmetic =', mpmath.nstr(Cth_aa_03, 15),
      ' run 9 published arith 2.51118750529894')
OUT['cert_Ctheta_aa_03'] = mpmath.nstr(Cth_aa_03, 15)

# ------------------------------------------------------------ the grid
EPS_GRID = [mpf(x) for x in ('0.30', '0.24', '0.19', '0.15', '0.12', '0.10')]

a_z = RHO[1]
b_z = 1 - a_z                      # = conj(rho_1)
ua_z, ub_z = u_of(a_z), u_of(b_z)
b_th = 1 - a_th                    # 1/4 - i gamma_1
ua_t, ub_t = u_of(a_th), u_of(b_th)
print('zeta labels: a =', mpmath.nstr(a_z, 15), ' b =', mpmath.nstr(b_z, 15))
print('theta labels: a =', mpmath.nstr(a_th, 15), ' b =', mpmath.nstr(b_th, 15))
OUT['labels'] = {'zeta_a': mpmath.nstr(a_z, 20), 'zeta_b': mpmath.nstr(b_z, 20),
                 'theta_a': mpmath.nstr(a_th, 20), 'theta_b': mpmath.nstr(b_th, 20)}

rows = []
for eps in EPS_GRID:
    # zeta pairs (arithmetic leg)
    Caa_z, _ = C_eps(eps, a_z, a_z, Lcache)
    Cab_z, _ = C_eps(eps, a_z, b_z, Lcache)
    Cbb_z, _ = C_eps(eps, b_z, b_z, Lcache)
    Dodd_z = Caa_z - 2 * Cab_z + Cbb_z
    Deven_z = Caa_z + 2 * Cab_z + Cbb_z
    Sodd_z = spec_odd_zeta(eps, ua_z, ub_z)
    # theta pairs (arithmetic leg)
    Caa_t, _ = C_theta(eps, a_th, a_th)
    Cab_t, _ = C_theta(eps, a_th, b_th)
    Cbb_t, _ = C_theta(eps, b_th, b_th)
    Dodd_t = Caa_t - 2 * Cab_t + Cbb_t
    Deven_t = Caa_t + 2 * Cab_t + Cbb_t
    Sodd_t = spec_odd_theta(eps, ua_t, ub_t)
    tail = spectral_tail_bound(eps)
    row = {'eps': str(eps),
           'zeta': {'Caa': mpmath.nstr(Caa_z, 15), 'Cab': mpmath.nstr(Cab_z, 15),
                    'Cbb': mpmath.nstr(Cbb_z, 15),
                    'Dodd_arith': mpmath.nstr(Dodd_z, 15),
                    'Deven_arith': mpmath.nstr(Deven_z, 15),
                    'Dodd_spectral': mpmath.nstr(Sodd_z, 15)},
           'theta': {'Caa': mpmath.nstr(Caa_t, 15), 'Cab': mpmath.nstr(Cab_t, 15),
                     'Cbb': mpmath.nstr(Cbb_t, 15),
                     'Dodd_arith': mpmath.nstr(Dodd_t, 15),
                     'Deven_arith': mpmath.nstr(Deven_t, 15),
                     'Dodd_spectral': mpmath.nstr(Sodd_t, 15)},
           'spectral_tail_bound_x3': mpmath.nstr(3 * tail, 4)}
    rows.append(row)
    print('eps', eps,
          '| Z odd arith', row['zeta']['Dodd_arith'],
          ' spec', row['zeta']['Dodd_spectral'],
          ' cross', row['zeta']['Cab'],
          '| TH odd arith', row['theta']['Dodd_arith'],
          ' spec', row['theta']['Dodd_spectral'],
          ' cross', row['theta']['Cab'], flush=True)
OUT['grid'] = rows

# ------------------------------------------------------------ extrapolation
import numpy as np


def fit_intercepts(eps_list, vals):
    out = {}
    xs = np.array([float(e) for e in eps_list])
    ys = np.array([float(v) for v in vals])
    for deg in (1, 2, 3):
        if len(xs) > deg:
            cf = np.polyfit(xs, ys, deg)
            out[f'deg{deg}'] = float(np.polyval(cf, 0.0))
    return out


def parse_real(s):
    return complex(s.replace(' ', '').replace('j', 'j')).real


for case in ('zeta', 'theta'):
    for qty in ('Dodd_arith', 'Dodd_spectral', 'Deven_arith', 'Cab'):
        vals = [parse_real(r[case][qty]) / (1 if qty == 'Cab' else 2) for r in rows]
        # G-normalisation: diagonals /2; cross reported as G too (C/2)
        if qty == 'Cab':
            vals = [parse_real(r[case][qty]) / 2 for r in rows]
        f_full = fit_intercepts(EPS_GRID, vals)
        sub = [(e, v) for e, v in zip(EPS_GRID, vals) if e <= mpf('0.15')]
        f_sub = fit_intercepts([e for e, _ in sub], [v for _, v in sub])
        allfits = dict(f_full)
        allfits.update({'sub_' + k: v for k, v in f_sub.items()})
        centrals = [v for k, v in allfits.items() if 'deg2' in k or 'deg3' in k]
        G0 = float(np.median(centrals)) if centrals else float('nan')
        err = (max(allfits.values()) - min(allfits.values())) / 2 if allfits else float('nan')
        key = f'G0_{case}_{qty}'
        OUT[key] = {'G0': G0, 'fit_spread_half': err, 'fits': allfits}
        print(key, '= ', G0, ' fit half-spread', err)

OUT['runtime_seconds'] = time.time() - T_START
with open('/home/hatch/workspace/gppdiscovery3/research/'
          'flint_odd_sector_2026_11_14_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('done', OUT['runtime_seconds'])
