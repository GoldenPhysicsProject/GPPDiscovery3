#!/usr/bin/env python
"""Flint tenth run (2026-11-14): verification of Heart's registered
prediction for the J-odd concentration diagonal of her sesquilinear
Weil form (her filing: GPP-bridge GROK_RESEARCH_NOTES.md @ 27691cf,
item 6; PR #11 comment 6888360390).

HER DEFINITIONS (verbatim).  Mellin transform fhat; the Weil form's
spectral face W(f,f) = sum_rho m_rho fhat(rho) conj(fhat(1 - conj(rho)));
J f(x) = x^{-1} conj(f(1/x)), so (Jf)^(s) = conj(fhat(1 - conj(s))), J
antilinear.  Concentration f = g_s - J g_s at label s, unit amplitude.

DICTIONARY (derived in-script, checked numerically at finite eps).
With tau(s) = 1 - conj(s) and the unit-peak Abel kernel
Ahat^q(w) = eps^2 / (eps^2 - (w - q)^2):
  g_s has transform Ahat^s;  (J g_s)^(w) = conj(Ahat^s(tau(w)))
  = Ahat^{tau(s)}(w)  (since tau(w) - s = -conj(w - tau(s)) up to the
  square), hence
  fhat(w) = Ahat^s(w) - Ahat^{tau(s)}(w),
  fhat(tau(w)) = -conj(fhat(w))   identically,
  W_eps(f, f) = - sum_w m_w |fhat(w)|^2   (manifestly <= 0).
The direct sum sum_w fhat(w) conj(fhat(tau(w))) is ALSO evaluated, as an
independent check of the identity.

REGISTERED PREDICTION (hers, before this run): W_eps -> -2 at
s_theta = rho_1 + 1/4 (off-line zero of xi_theta, m = 1, a = 1), and
W_eps = 0 at rho_1 for zeta (on-line: tau(s) = s, so fhat = 0 identically).

WHY SPECTRAL LEG ONLY: W is sesquilinear; its arithmetic face would need
the corpus evaluator to accept an antiholomorphic second slot, which it
cannot -- Heart's item 5 obstruction, in evaluator form.  The spectral
leg is the falsifier evaluator's own certified leg (run 9: the two legs
agreed to 5e-7 at the theta diagonal).

Variants reported for the dictionary table (same sums, other
normalisations): peak-2 corpus kernels without the 1/2 (limit -8) and
with the corpus 1/2 (limit -4).  Her '-2 (G-normalised -1)' is the
unit-amplitude value and its half.

Epistemic status: NUMERICAL; modulo TailLemma for the identification of
the Abel-regularised limit with her idealised two-point concentration,
exactly as in her filing.
"""

import json
import sys
import time

sys.path.insert(0, '/home/hatch/workspace/zeros-platt')

import mpmath
from mpmath import mp, mpf, mpc
mp.dps = 50

from platt_sequential import zeros_in_range

T_START = time.time()
OUT = {}
GAM1 = mpf('14.13472514173469379045725198356247027078425711569924328783754417118910251371121454762')
RHO1 = mpc(mpf('0.5'), GAM1)


def tau(s):
    return 1 - s.conjugate()


def Ahat(q, w, eps):
    return eps ** 2 / (eps ** 2 - (w - q) ** 2)


def W_jodd(eps, s, theta, T=3000):
    """Returns (W via -sum|fhat|^2, W direct, tail bound)."""
    ts = tau(s)
    tot_id = mpc(0)
    tot_dir = mpc(0)
    gammas = zeros_in_range(0, T)
    for g in gammas:
        gm = mpf(g)
        if theta:
            ws = [mpc(mpf('0.25'), gm), mpc(mpf('0.25'), -gm),
                  mpc(mpf('0.75'), gm), mpc(mpf('0.75'), -gm)]
        else:
            ws = [mpc(mpf('0.5'), gm), mpc(mpf('0.5'), -gm)]
        for w in ws:
            fw = Ahat(s, w, eps) - Ahat(ts, w, eps)
            tot_id -= abs(fw) ** 2
            ftw = Ahat(s, tau(w), eps) - Ahat(ts, tau(w), eps)
            tot_dir += fw * ftw.conjugate()
    # tail majorant: |Ahat^q(w)| <= eps^2 / ((|Im w| - |Im q| - 1)^2), so
    # |fhat|^2 <= 2(|Ahat^s|^2 + |Ahat^ts|^2) <= 4 eps^4 / (t - g1 - 1)^4;
    # zero density (log t)/(2 pi) per unit height, both signs, x2 families
    # for theta (factor 4 total vs zeta's 2) -- use the generous factor 4.
    from mpmath import quad, log, inf, pi
    tail = quad(lambda t: 4 * (log(t) / (2 * pi)) * 4 * eps ** 4 / (t - 16) ** 4,
                [T, inf])
    return tot_id, tot_dir, tail


EPS_GRID = [mpf(x) for x in ('0.30', '0.24', '0.19', '0.15', '0.12', '0.10',
                             '0.08', '0.06')]
for case, s, theta in (('zeta_rho1', RHO1, False),
                       ('theta_stheta', RHO1 + mpf('0.25'), True)):
    rows = []
    for eps in EPS_GRID:
        wid, wdir, tail = W_jodd(eps, s, theta)
        rows.append({'eps': str(eps), 'W_identity': mpmath.nstr(wid, 15),
                     'W_direct': mpmath.nstr(wdir, 15),
                     'tail_bound': mpmath.nstr(tail, 3)})
        print(case, 'eps', eps, ' W =', rows[-1]['W_identity'],
              ' direct =', rows[-1]['W_direct'], ' tail <=', rows[-1]['tail_bound'],
              flush=True)
    OUT[case] = rows
    import numpy as np
    xs = np.array([float(r['eps']) for r in rows])
    ys = np.array([float(complex(r['W_identity'].strip('()').replace(' ', '')).real)
                   for r in rows])
    fits = {}
    for deg in (1, 2, 3):
        cf = np.polyfit(xs, ys, deg)
        fits[f'deg{deg}'] = float(np.polyval(cf, 0.0))
    centrals = [fits['deg2'], fits['deg3']]
    OUT[case + '_extrapolation'] = {
        'fits': fits, 'W0': float(np.median(centrals)),
        'fit_spread_half': (max(fits.values()) - min(fits.values())) / 2}
    print(case, 'extrapolated W0 =', OUT[case + '_extrapolation'], flush=True)

OUT['runtime_seconds'] = time.time() - T_START
with open('/home/hatch/workspace/gppdiscovery3/research/'
          'flint_heart_jodd_2026_11_14_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('done', OUT['runtime_seconds'])
