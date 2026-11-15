#!/usr/bin/env python
"""Flint eleventh run (2026-11-15): the J-Gram from Flint's own machinery.

CONTEXT.  Run 10's queue item 10 asked for all four entries of the
2x2 Gram matrix of (g_s, J g_s) in Heart's sesquilinear Weil pairing
    W(f, h) = sum_w m_w f^(w) * conj( h^(tau w) ),   tau w = 1 - conj(w),
at the two test labels: zeta's first zero rho_1 and F_theta's first
zero s_theta = rho_1 + 1/4 (theta = 1/4, F_theta(w) = xi(w-1/4) xi(w+1/4)).
Heart answered Q3 in her Rounds 2-3 (GPPDiscovery5 0ff0650, 950b1b3):
the entries are the absolutely convergent zero sums
    B(q, e) = sum_w m_w Ahat^q(w)^2,
    A(q, e) = sum_w m_w Ahat^q(w) Ahat^{tau q}(w),
    Ahat^q(w) = e^2 / (e^2 - (w - q)^2)     (unit peak),
with W(g, Jg) = B, W(g, g) = A, W(Jg, Jg) = A, W(Jg, g) = conj(B),
D := 2 Re A - 2 Re B (the J-odd diagonal), det = A^2 - |B|^2.
Her J2 evaluates them by residues of L = F'/F at q+-e, tau q+-e.

THIS SCRIPT computes the entries DEFINITION-SIDE: direct zero sums
over the staged Platt ordinates (|gamma| <= 3000, multiplicities 1,
the convention of Heart's own 400-zero sums), at dps 50, for
e in {0.12, 0.06, 0.03}; plus an INDEPENDENT implementation of her
J2 residue formulas (L from mpmath zeta and digamma) at the same
epsilons, so the smaller epsilon 0.03 -- where Heart filed no
values -- is covered by two legs of mine.  Filed comparison targets
(Heart Round 3, GPPDiscovery5 950b1b3, jgram_residue_results.json;
Round 4 statement sheet J4):
  rho_1, e=0.06:  A = B = 1.00000000715947, D = 0, det = 0.
  s_theta, e=0.06: A = -0.0292207650958216,
                   B = 1.00021347748317 - 1.83693300282492e-9 i,
                   D = -2.05886848515798, det = -0.999573147426184.
  s_theta, e=0.12: A = -0.1222408606753, B = 1.003735944705 - 2.94e-8 i,
                   D = -2.25195361076, det = -0.9925430186734.

Also computed: B(tau s_theta, e) over F_theta zeros (the shadow
argument -- run 9's completed diagonal G_e(s_theta) is this number;
see flint_mevaluator_compare_2026_11_15.py), and direct-sum checks
of W(Jg, Jg) = A and W(Jg, g) = conj(B) from (Jg)^(w) := conj(Ahat^q(tau w)).

Tail control: far zeros |gamma| > 3000 contribute, by the joint
TailLemma T0/T2 majorant (density (1/2pi) log(t/2pi) + 1),
<= (16/9) e^4 * S where S is the majorant integral evaluated
in-script; printed per row.  Zeros are used only as the summation
set of the DEFINITION of A and B (and as test labels), never as
input to any construction.

Epistemic status: NUMERICAL (dps 50 zero sums + residue leg).
"""

import json
import sys
import time

sys.path.insert(0, '/home/hatch/workspace/zeros-platt')

import mpmath
from mpmath import mp, mpf, mpc, pi, log, digamma, polygamma, quad

mp.dps = 50
T_START = time.time()
OUT = {'labels': {}, 'tail_majorants': {}}

from platt_sequential import zeros_in_range, zeros_by_index

GAM = sorted(zeros_in_range(0, 3000))
GAM1 = mpf(dict(zeros_by_index(1, 1))[1])
print('ordinates staged:', len(GAM), ' gamma_1 =', mpmath.nstr(GAM1, 30))
OUT['n_ordinates'] = len(GAM)
OUT['gamma_1'] = mpmath.nstr(GAM1, 40)

TH = mpf('0.25')
RHO1 = mpc(mpf('0.5'), GAM1)
S_TH = mpc(mpf('0.75'), GAM1)          # rho_1 + 1/4, zero of F_theta
TAU_S_TH = mpc(mpf('0.25'), GAM1)      # tau(s_theta) = rho_1 - 1/4


def Ahat(q, w, e):
    return e**2 / (e**2 - (w - q)**2)


def tau(w):
    return 1 - w.conjugate()


# zero sets as (w, m) lists
Z_XI = []
for g in GAM:
    Z_XI.append(mpc(mpf('0.5'), g))
    Z_XI.append(mpc(mpf('0.5'), -g))
Z_FTH = []
for g in GAM:
    for re in (mpf('0.25'), mpf('0.75')):
        Z_FTH.append(mpc(re, g))
        Z_FTH.append(mpc(re, -g))


def gram_sums(Z, q, e):
    """A, B, W(Jg,Jg), W(Jg,g) by direct definition-side sums."""
    tq = tau(q)
    A = mpc(0)
    B = mpc(0)
    JJ = mpc(0)
    Jg_ = mpc(0)
    for w in Z:
        a_q = Ahat(q, w, e)
        a_tq = Ahat(tq, w, e)
        tw = tau(w)
        a_q_tw = Ahat(q, tw, e)          # for (Jg)^(w) = conj(Ahat^q(tau w))
        A += a_q * a_tq
        B += a_q * a_q
        jg_w = a_q_tw.conjugate()        # (Jg)^(w)
        # (Jg)^(tau w) = conj(Ahat^q(tau tau w)) = conj(Ahat^q(w))
        JJ += jg_w * (a_q.conjugate()).conjugate()
        # W(Jg, g): sum (Jg)^(w) * conj( g^(tau w) ), g^(tau w) = Ahat^q(tau w)
        Jg_ += jg_w * a_q_tw.conjugate()
    return A, B, JJ, Jg_


def tail_majorant(e, gam_q, n_series=2):
    """Majorant for the |gamma| > 3000 tail of a Gram-type sum:
    (16/9) e^4 * int_{3000}^inf [(1/2pi) log(t/2pi) + 1] (t - gam_q)^-4 dt
    times n_series (the number of zero series: 2 for xi (+-gamma),
    4 for F_theta)."""
    f = lambda t: ((log(t / (2 * pi)) / (2 * pi)) + 1) / (t - gam_q) ** 4
    return n_series * (mpf(16) / 9) * e**4 * quad(f, [3000, mpmath.inf])


# ------------------------------------------------------------- J2 leg
def L_xi(w):
    z = mp.zeta(w)
    zp = mp.zeta(w, derivative=1)
    return (zp / z + 1 / w + 1 / (w - 1)
            - log(pi) / 2 + digamma(w / 2) / 2)


def Lp_xi(w):
    z = mp.zeta(w)
    zp = mp.zeta(w, derivative=1)
    zpp = mp.zeta(w, derivative=2)
    return (zpp / z - (zp / z) ** 2 - 1 / w**2 - 1 / (w - 1) ** 2
            + polygamma(1, w / 2) / 4)


def L_th(w):
    return L_xi(w - TH) + L_xi(w + TH)


def Lp_th(w):
    return Lp_xi(w - TH) + Lp_xi(w + TH)


def B_residue(L, Lp, q, e):
    return -((e**2 / 4) * (Lp(q + e) + Lp(q - e))
             + (e / 4) * (L(q - e) - L(q + e)))


def A_residue(L, q, e):
    tq = tau(q)
    if tq == q:
        return None  # A = B identically
    return -((-(e / 2) * Ahat(tq, q + e, e) * L(q + e)
              + (e / 2) * Ahat(tq, q - e, e) * L(q - e)
              - (e / 2) * Ahat(q, tq + e, e) * L(tq + e)
              + (e / 2) * Ahat(q, tq - e, e) * L(tq - e)))


FILED = {
    ('xi', '0.06'): {'A': '1.00000000715947', 'B': '1.00000000715947',
                     'D': '0', 'det': '0'},
    ('fth', '0.06'): {'A': '-0.0292207650958216',
                      'B': '(1.00021347748317 - 1.83693300282492e-9j)',
                      'D': '-2.05886848515798', 'det': '-0.999573147426184'},
    ('fth', '0.12'): {'A': '-0.1222408606753',
                      'B': '(1.003735944705 - 2.937274374611e-8j)',
                      'D': '-2.25195361076', 'det': '-0.9925430186734'},
}

CASES = [
    ('xi', Z_XI, RHO1, L_xi, Lp_xi, 2),
    ('fth', Z_FTH, S_TH, L_th, Lp_th, 4),
]

for name, Z, q, L, Lp, nser in CASES:
    for e in (mpf('0.12'), mpf('0.06'), mpf('0.03')):
        A, B, JJ, Jg_ = gram_sums(Z, q, e)
        D = 2 * A.real - 2 * B.real
        det = A * A - B * B.conjugate()
        Br = B_residue(L, Lp, q, e)
        Ar = A_residue(L, q, e)
        if Ar is None:
            Ar = Br
        tm = tail_majorant(e, abs(q.imag), nser)
        row = {
            'eps': str(e),
            'A_sum': mpmath.nstr(A, 25), 'B_sum': mpmath.nstr(B, 25),
            'WJJ_sum': mpmath.nstr(JJ, 20), 'WJg_sum': mpmath.nstr(Jg_, 20),
            'D_from_sums': mpmath.nstr(D, 25),
            'det_from_sums': mpmath.nstr(det, 25),
            'D_over_B': mpmath.nstr(D / B.real, 20),
            'A_residue_mine': mpmath.nstr(Ar, 25),
            'B_residue_mine': mpmath.nstr(Br, 25),
            'sum_minus_res_A': mpmath.nstr(abs(A - Ar), 4),
            'sum_minus_res_B': mpmath.nstr(abs(B - Br), 4),
            'Im_A': mpmath.nstr(abs(A.imag), 4),
            'WJJ_minus_A': mpmath.nstr(abs(JJ - A), 4),
            'WJg_minus_conjB': mpmath.nstr(abs(Jg_ - B.conjugate()), 4),
            'tail_majorant': mpmath.nstr(tm, 4),
        }
        key = (name, str(e))
        if key in FILED:
            f = FILED[key]
            row['filed'] = f
            fb = complex(f['B'].strip('()').replace(' ', ''))
            row['sum_minus_filed_A'] = mpmath.nstr(abs(A - mpc(f['A'])), 4)
            row['sum_minus_filed_B'] = mpmath.nstr(
                abs(B - mpc(fb.real, fb.imag)), 4)
            row['sum_minus_filed_D'] = mpmath.nstr(abs(D - mpf(f['D'])), 4)
            row['sum_minus_filed_det'] = mpmath.nstr(
                abs(det - mpc(f['det'])), 4)
        OUT['labels'].setdefault(name, []).append(row)
        print(name, 'eps', e, 'A', row['A_sum'][:28], ' B', row['B_sum'][:28],
              ' D', row['D_from_sums'][:20], ' det', row['det_from_sums'][:22],
              ' |sum-res| B', row['sum_minus_res_B'],
              (' |sum-filed| B ' + row['sum_minus_filed_B'])
              if 'sum_minus_filed_B' in row else '')

# B at the shadow label of F_theta (task 3 ground truth, definition side)
for e in (mpf('0.12'), mpf('0.06')):
    Bs = mpc(0)
    for w in Z_FTH:
        a = Ahat(TAU_S_TH, w, e)
        Bs += a * a
    Br = B_residue(L_th, Lp_th, TAU_S_TH, e)
    OUT.setdefault('B_at_tau_s_theta', []).append({
        'eps': str(e), 'B_sum': mpmath.nstr(Bs, 25),
        'B_residue_mine': mpmath.nstr(Br, 25),
        'diff': mpmath.nstr(abs(Bs - Br), 4)})
    print('B(tau s_theta) eps', e, mpmath.nstr(Bs, 25),
          ' residue', mpmath.nstr(Br, 25), ' diff', mpmath.nstr(abs(Bs - Br), 4))

OUT['runtime_seconds'] = time.time() - T_START
with open('/home/hatch/workspace/gppdiscovery3/research/'
          'flint_jgram_2026_11_15_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('done', OUT['runtime_seconds'])
