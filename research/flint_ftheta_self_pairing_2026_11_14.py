#!/usr/bin/env python
"""Flint ninth run (2026-11-14): the F_theta SELF-pairing concentration.

QUESTION.  The eighth run's extracted-state norm formula says the state
extracted at a zero rho of the paired function has completed Gram norm
    ||[Psi_rho]||^2 = m_rho * lim_{eps->0} eps^2/(eps^2 - (2 sigma_rho - 1)^2),
i.e. a null class iff the zero is off the line.  That formula was assembled
from the residue of the completed pairing P(a) = -(xi'/xi)(1/2 + a)
(run 8, C1) times the FREE diagonal Cesaro factor (C3).  It was never
tested against the completed pairing of a function that actually HAS an
off-line zero.  The falsifier arm tested (a) zeta's completed pairing at
zeta's on-line zeros (Gram -> I, SURVIVES) and (b) zeta's completed
pairing at the SHIFTED labels rho +- 1/4 (-> 0, because those labels are
not zeros of zeta: their kernel centres miss zeta's zero set).  The
discriminating case is (c): F_theta's OWN completed pairing at F_theta's
OWN zero.  This script computes exactly that.

OBJECT.  xi_theta(s) = xi(s - 1/4) xi(s + 1/4); its zeros are
rho' +- 1/4 for zeta zeros rho'.  First theta zero: s_th = rho_1 + 1/4
(sigma = 3/4).  By Hadamard for xi_theta (entire, order 1, FE
xi_theta(s) = xi_theta(1 - s); THEOREM, classical),
    W_theta[F] := sum_{rho_theta} F(rho_theta)
               = W[F(. - 1/4)] + W[F(. + 1/4)],
with W the Guinand-Weil functional of cesaro_gram_falsifier_2026_11_14.py.
A spectral weight centred at u becomes centre u +- 1/4 under the shifts,
so with u_th = 1 - conj(s_th) = 1/4 + i gamma_1:
    C^theta_eps = Cmach(eps; u_th - 1/4) + Cmach(eps; u_th + 1/4),
where Cmach(eps; u') is the falsifier's arithmetic assembly
(pole - PRIME by continuation + ARCH) for the weight with centre u'.

PREDICTION (registered in rh_lab/catalog.json before the run): the
spectral sum's term at the SHADOW zero rho_1 - 1/4 = 1/4 + i gamma_1 sits
exactly at the centre u_th, contributing full weight 2; all other terms
vanish as eps -> 0.  Hence C^theta_eps -> 2 (corpus), G^theta = C^theta/2
-> 1 -- the same value as zeta on-line -- and NOT the free-selector
value 2 eps^2/(eps^2 - 1/4) -> 0 that the eighth run's formula predicts.
If the extrapolation instead lands on 0 within bars, the prediction is
wrong and the eighth run's formula survives.

EVALUATOR.  The falsifier script's definitions are exec'd byte-for-byte
(RUN section stripped); before any theta number, the evaluator is
re-certified on three published values of that script: C(0.3; rho_1,
rho_1) arithmetic = 2.0000089211, its spectral sum, and the shifted
zeta diagonal C(0.3; rho_1 + 1/4, rho_1 + 1/4) from its results JSON.

Epistemic status: NUMERICAL.  The spectral-face identity for xi_theta
is a theorem; the run certifies the limit and the evaluator.
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
W_val = NS['W_val']; u_of = NS['u_of']
_quad = NS['quad']; _log = NS['log']

def spectral_tail_bound(eps, T=3000):
    """The falsifier script's spectral tail majorant (its RUN section,
    quoted verbatim): bound for the |gamma| > T tail of the spectral sum."""
    return _quad(lambda t: 2 * eps ** 4 * (_log(mpf(t)) + 13) / (mpf(t) - 30) ** 4,
                 [T, mpmath.inf])

T_START = time.time()
OUT = {'evaluator': 'cesaro_gram_falsifier_2026_11_14.py definitions, exec\'d verbatim'}

Lcache = {}

# ------------------------------------------------------------ certification
print('== Certification of the reused evaluator ==')
C11_ar, _ = C_eps(mpf('0.3'), RHO[1], RHO[1], Lcache)
print('C(0.3; rho1, rho1) arithmetic =', mpmath.nstr(C11_ar, 18),
      ' published 2.0000089211')
OUT['cert_C11_arith'] = mpmath.nstr(C11_ar, 18)
S11 = spectral_sum(mpf('0.3'), u_of(RHO[1]), u_of(RHO[1]))
print('C(0.3; rho1, rho1) spectral   =', mpmath.nstr(S11, 18),
      ' tail bound', mpmath.nstr(spectral_tail_bound(mpf('0.3')), 3))
OUT['cert_C11_spectral'] = mpmath.nstr(S11, 18)

pub = json.load(open('/home/hatch/workspace/gppdiscovery3/research/'
                     'cesaro_gram_falsifier_2026_11_14_results.json'))
pub_shift = pub['control_completed_shift']['k=1 shift=+0.25']['0.3']
s_shift = RHO[1] + mpf('0.25')
Csh_ar, _ = C_eps(mpf('0.3'), s_shift, s_shift, Lcache)
print('C(0.3; rho1+1/4) zeta-pairing arithmetic =', mpmath.nstr(Csh_ar, 15),
      ' published', pub_shift)
OUT['cert_shift_zeta'] = {'arith': mpmath.nstr(Csh_ar, 15), 'published': pub_shift}

# ------------------------------------------------------------ theta run
def s_fake(u_prime):
    """Fake label whose amplitude centre u_of(s) equals u_prime."""
    return 1 - u_prime.conjugate()

def C_theta_arith(eps, u_th):
    """W_theta[W_{u_th}] via the two shifted zeta assemblies."""
    tot = mpc(0)
    parts = {}
    for sgn, up in (('-1/4', u_th - mpf('0.25')), ('+1/4', u_th + mpf('0.25'))):
        v, det = C_eps(eps, s_fake(up), s_fake(up), Lcache)
        parts[sgn] = mpmath.nstr(v, 15)
        tot += v
    return tot, parts

def C_theta_spectral(eps, u_th, T=3000):
    from platt_sequential import zeros_in_range
    S = mpc(0)
    main = None
    for g in zeros_in_range(0, T):
        for re in (mpf('0.25'), mpf('0.75')):
            for sgn in (1, -1):
                w = mpc(re, sgn * g)
                val = W_val(w, eps, u_th, u_th)
                S += val
                if re == mpf('0.25') and sgn == 1 and g == GAM[1]:
                    main = val
    return S, main

s_th = RHO[1] + mpf('0.25')          # 3/4 + i gamma_1, a zero of xi_theta
u_th = u_of(s_th)                    # 1/4 + i gamma_1
print('theta label s_th =', mpmath.nstr(s_th, 20), ' centre u_th =',
      mpmath.nstr(u_th, 20))
OUT['s_theta'] = mpmath.nstr(s_th, 20)

EPS_GRID = [mpf(x) for x in ('0.30', '0.26', '0.225', '0.195', '0.17',
                             '0.15', '0.13', '0.11', '0.09', '0.07')]
rows = []
for eps in EPS_GRID:
    Ca, parts = C_theta_arith(eps, u_th)
    Cs, main = C_theta_spectral(eps, u_th)
    tail = spectral_tail_bound(eps)
    freefac = 2 * eps ** 2 / (eps ** 2 - mpf('0.25'))
    row = {'eps': str(eps),
           'C_theta_arith': mpmath.nstr(Ca, 15),
           'C_theta_spectral': mpmath.nstr(Cs, 15),
           'legs_diff': mpmath.nstr(abs(Ca - Cs), 4),
           'spectral_tail_bound': mpmath.nstr(tail, 4),
           'shadow_main_term': mpmath.nstr(main, 15),
           'rest_of_spectral': mpmath.nstr(Cs - main, 15),
           'shift_parts_arith': parts,
           'run8_formula_prediction_C': mpmath.nstr(freefac, 12)}
    rows.append(row)
    print('eps', eps, 'arith', row['C_theta_arith'], ' spec',
          row['C_theta_spectral'], ' diff', row['legs_diff'],
          ' main', row['shadow_main_term'], ' rest', row['rest_of_spectral'],
          ' run8-pred', row['run8_formula_prediction_C'])
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

for leg in ('C_theta_arith', 'C_theta_spectral'):
    vals = []
    for r in rows:
        z = complex(r[leg].replace(' ', '').replace('j', 'j'))
        vals.append(z.real / 2)
    f_full = fit_intercepts(EPS_GRID, vals)
    sub = [(e, v) for e, v in zip(EPS_GRID, vals) if e <= mpf('0.17')]
    f_sub = fit_intercepts([e for e, _ in sub], [v for _, v in sub])
    allfits = dict(f_full)
    allfits.update({'sub_' + k: v for k, v in f_sub.items()})
    centrals = [v for k, v in allfits.items() if 'deg2' in k or 'deg3' in k]
    G0 = float(np.median(centrals))
    err = (max(allfits.values()) - min(allfits.values())) / 2
    OUT[f'G0_theta_{leg}'] = {'G0': G0, 'fit_spread_half': err,
                              'fits': allfits}
    print(leg, 'G0^theta =', G0, ' fit half-spread', err, allfits)

OUT['runtime_seconds'] = time.time() - T_START
with open('/home/hatch/workspace/gppdiscovery3/research/'
          'flint_ftheta_self_pairing_2026_11_14_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('done', OUT['runtime_seconds'])
