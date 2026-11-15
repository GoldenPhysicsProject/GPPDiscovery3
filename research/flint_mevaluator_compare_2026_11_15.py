#!/usr/bin/env python
"""Flint eleventh run (2026-11-15): run 9's M evaluator vs Heart's J2.

QUESTION (Heart Round 3, J5 follow-up; Turbo's audit 6889087427 item
(3) leaves it as the open comparison).  Run 9 evaluated F_theta's
completed pairing at its own zero s_theta by the falsifier assembly:
    C^theta_e = W[F(. - 1/4)] + W[F(. + 1/4)]   (Hadamard shift),
each W[.] assembled arithmetic-side (pole - PRIME by continuation +
ARCH) for the corpus weight W_{kk,e}(w) = (1/2) A_e(s, w)^2 with
A_e(s, w) = 2 e^2/(e^2 - (w - u)^2), u = 1 - conj(s) (peak 2,
SHADOW centre).  Heart's J-Gram B entry is the unit-peak square sum
at a centre q:  B(q, e) = sum_w m_w (Ahat^q(w))^2.

ALGEBRA (checked in-script, pointwise): A_e(s, w) = 2 Ahat^{tau s}(w)
identically, so (1/2) A_e(s, .)^2 = 2 (Ahat^{tau s})^2 and therefore,
on spectral faces, for either F = xi or F = F_theta,
    G_e(s) := C_e(s)/2  =  B(tau s, e)
EXACTLY -- run 9's "M evaluator" is Heart's B read at the SHADOW
argument.  What is NOT algebra is whether run 9's arithmetic
ASSEMBLY (continuation in e from the Weil-class region e > 1/2,
per-face: PRIME via L(1/2 -+ a), ARCH by quadrature at off-line
centres) computes the same number as the zero sum / J2 residue form
in the concentrating regime e < 1/2, where the presentations could
in principle separate (J5: no convergent prime sum there).

THIS SCRIPT computes the assembly at e = 0.12 and e = 0.06 for the
F_theta diagonal at s_theta (two shifted C_eps calls, falsifier
definitions exec'd verbatim, exactly as run 9) and for the zeta
on-line diagonal at rho_1, and compares against the definition-side
zero sums and J2 values from flint_jgram_2026_11_15.py:
  B(tau s_theta, 0.12) = 1.0037359447046086 + 2.93727437461e-8 i
  B(tau s_theta, 0.06) = 1.0002134774831662 + 1.83693300282e-9 i
  B(rho_1, 0.12)       = 1.00000011450638703 (real; tau rho_1 = rho_1)

Epistemic status: NUMERICAL.  The assembly leg carries the
falsifier's own error budget (Arb L values prec 250; arch tail
T = 400 analytic bound; run 9 filed an e-dependent arch-quadrature
error at off-line centres as its leg-gap caveat).
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
C_eps = NS['C_eps']

T_START = time.time()
OUT = {}

# ------------------------------------------------ algebra identity
e0 = mpf('0.06')
s0 = mpc(mpf('0.75'), GAM[1])
u0 = 1 - s0.conjugate()
worst = mpf(0)
for wtest in (mpc(mpf('0.3'), mpf('40.2')), mpc(mpf('0.9'), mpf('-7.7')),
              mpc(mpf('0.25'), GAM[1])):
    A_e = 2 * e0**2 / (e0**2 - (wtest - u0) ** 2)
    Ahat = e0**2 / (e0**2 - (wtest - u0) ** 2)
    worst = max(worst, abs(A_e / 2 - Ahat))
OUT['algebra_half_Aeps_eq_Ahat'] = mpmath.nstr(worst, 4)
print('algebra: max |(1/2) A_e - Ahat^{tau s}| at sample points =',
      OUT['algebra_half_Aeps_eq_Ahat'])

Lcache = {}


def s_fake(u_prime):
    return 1 - u_prime.conjugate()


def C_theta_arith(eps, u_th):
    tot = mpc(0)
    parts = {}
    for sgn, up in (('-1/4', u_th - mpf('0.25')), ('+1/4', u_th + mpf('0.25'))):
        v, det = C_eps(eps, s_fake(up), s_fake(up), Lcache)
        parts[sgn] = mpmath.nstr(v, 18)
        tot += v
    return tot, parts


s_th = RHO[1] + mpf('0.25')
u_th = NS['u_of'](s_th)
print('s_theta =', mpmath.nstr(s_th, 20), ' u_th =', mpmath.nstr(u_th, 20))

TARGETS_TH = {
    '0.12': mpc(mpf('1.003735944704608615663041'),
                mpf('2.937274374605970542614665e-8')),
    '0.06': mpc(mpf('1.000213477483166237323675'),
                mpf('1.836933002821866563821482e-9')),
}
rows = []
for eps in (mpf('0.12'), mpf('0.06')):
    Ca, parts = C_theta_arith(eps, u_th)
    G = Ca / 2
    tgt = TARGETS_TH[str(eps)]
    row = {'eps': str(eps), 'assembly_G': mpmath.nstr(G, 20),
           'assembly_parts': parts,
           'B_tau_stheta_zerosum': mpmath.nstr(tgt, 20),
           'assembly_minus_zerosum': mpmath.nstr(abs(G - tgt), 4),
           'rel': mpmath.nstr(abs(G - tgt) / abs(tgt), 4)}
    rows.append(row)
    print('theta eps', eps, 'assembly G =', row['assembly_G'],
          ' |assembly - B(tau s_theta)| =', row['assembly_minus_zerosum'])
OUT['theta_rows'] = rows

# zeta on-line diagonal: assembly C_eps(rho_1, rho_1)/2 vs B(rho_1, e)
Cz, _ = C_eps(mpf('0.12'), RHO[1], RHO[1], Lcache)
Gz = Cz / 2
tgtz = mpf('1.00000011450638703165198')
OUT['zeta_rho1_eps012'] = {
    'assembly_G': mpmath.nstr(Gz, 20),
    'B_rho1_zerosum': mpmath.nstr(tgtz, 20),
    'assembly_minus_zerosum': mpmath.nstr(abs(Gz - tgtz), 4)}
print('zeta rho_1 eps 0.12 assembly G =', OUT['zeta_rho1_eps012']['assembly_G'],
      ' |diff| =', OUT['zeta_rho1_eps012']['assembly_minus_zerosum'])

OUT['runtime_seconds'] = time.time() - T_START
with open('/home/hatch/workspace/gppdiscovery3/research/'
          'flint_mevaluator_compare_2026_11_15_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('done', OUT['runtime_seconds'])
