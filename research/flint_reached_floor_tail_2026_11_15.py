#!/usr/bin/env python
"""Run 14 supplement to run 13 (flint_reached_floor_2026_11_15):
derivation of the (H-Env) tail input used in the assembled bound
of run 13, section 2 (P5).

Run 13's assembled bound consumed the tail input
    T := sum_{j>812} y_j^2 / A_jj  <=  4.8991e-12
with no derivation in the committed artifacts (Turbo audit, run 14).
This script derives it from the two named hypotheses, using only
constants already in run 13's committed results JSON:

  (H-Env)  |y_j| <= C_y * omega_j^-3   for all j >= 73,
           C_y = 3.3293405873747175e-3 (measured max on 73..812;
           JSON reached_profile.y_v0_envelope_max_abs_yj_times_omj3)
  (H-A)    A_jj >= c_A * log(omega_j) / omega_j^4   for all j,
           c_A = 1.7372355768219303 (measured profile minimum on Z;
           JSON proof_quantities.A_diag_profile_om4_over_logom_min)

with omega_j = 4*pi*j.  Hence for every j > 812,
    y_j^2 / A_jj <= (C_y^2 / c_A) * 1 / (omega_j^2 * log omega_j),
and T <= (C_y^2 / c_A) * S with
    S := sum_{j>=813} 1 / (omega_j^2 * log omega_j).

S is evaluated as an explicit float64 partial sum to N = 200000
plus a rigorous integral bracket for the remainder: f(x) =
1/((4*pi*x)^2 * log(4*pi*x)) is decreasing for x >= 1, so
    int_{N+1}^inf f  <=  sum_{j>N} f(j)  <=  f(N+1) + int_{N+1}^inf f,
the integrals evaluated with mpmath at 30 digits.  The bracket is
widened by 1e-9 relative to cover the float64 partial-sum error
(~1e-13 relative at these magnitudes).

Outputs the tail bracket, the bound actually used (bracket upper
end), and the reassembled v0^T L v0 bound with its room vs beta0.
Conditional on (H-Env) and (H-A), exactly as run 13 consumes them.
No zero data.  Epistemic status: NUMERICAL derivation of a
conditional analytic bound.
"""

import json
import math

import mpmath as mp

mp.mp.dps = 30

SRC = ('/home/hatch/workspace/gppdiscovery3/research/'
       'flint_reached_floor_2026_11_15_results.json')
R13 = json.load(open(SRC))

C_Y = R13['reached_profile']['y_v0_envelope_max_abs_yj_times_omj3']
C_A = R13['proof_quantities']['A_diag_profile_om4_over_logom_min']
RHO_Z = R13['proof_quantities']['rho_Z_Dnormalised_offdiag_opnorm']
SIGMA_V0 = R13['thresholds']['SigmaV0_norm_vtildeT_v0_sq']
RQ_D = R13['profile_route']['RQ_D_sum_y2_over_Ajj_per_unit_norm']
BETA0 = R13['thresholds']['beta0']
BETA1 = R13['thresholds']['beta1']
FILED_TAIL = 4.8991e-12          # run 13 note, section 2 (P5), as filed

N = 200000


def f_float(j):
    w = 4.0 * math.pi * j
    return 1.0 / (w * w * math.log(w))


partial = 0.0
for j in range(813, N + 1):
    partial += f_float(j)


def f_mp(x):
    w = 4 * mp.pi * x
    return 1 / (w ** 2 * mp.log(w))


int_tail = mp.quad(f_mp, [N + 1, mp.inf])     # int_{N+1}^inf f
f_next = f_mp(N + 1)
S_lo = partial + float(int_tail)
S_hi = partial + float(f_next) + float(int_tail)
# widen for float64 partial-sum error
S_lo *= (1 - 1e-9)
S_hi *= (1 + 1e-9)

factor = C_Y ** 2 / C_A
tail_lo, tail_hi = factor * S_lo, factor * S_hi
finite = RQ_D * SIGMA_V0                       # sum_{j<=812} y_j^2/A_jj
total_hi = finite + tail_hi
assembled = total_hi / (1.0 - RHO_Z)

OUT = {
    'label': 'NUMERICAL derivation of the conditional (H-Env) tail '
             'bound for run 13 P5; hypotheses (H-Env), (H-A)',
    'inputs_from_run13_json': {
        'C_y_envelope': C_Y, 'c_A_profile_min': C_A,
        'rho_Z': RHO_Z, 'SigmaV0': SIGMA_V0, 'RQ_D': RQ_D,
        'beta0': BETA0, 'beta1': BETA1},
    'S_sum_bracket': {'partial_sum_813_to_N': partial, 'N': N,
                      'integral_bracket_remainder': float(int_tail),
                      'S_lo': S_lo, 'S_hi': S_hi},
    'factor_Cy2_over_cA': factor,
    'tail_bracket': {'lo': tail_lo, 'hi': tail_hi},
    'filed_tail_input_run13': FILED_TAIL,
    'filed_tail_inside_derived_bracket': bool(
        tail_lo <= FILED_TAIL <= tail_hi),
    'tail_bound_used': tail_hi,
    'finite_part_sum_y2_over_Ajj_j_le_812': finite,
    'total_bound': total_hi,
    'assembled_v0Lv0_bound': assembled,
    'room_vs_beta0': BETA0 / assembled,
    'room_vs_beta1': BETA1 / assembled,
}

with open('/home/hatch/workspace/gppdiscovery3/research/'
          'flint_reached_floor_tail_2026_11_15_results.json', 'w') as fh:
    json.dump(OUT, fh, indent=1)
print(json.dumps(OUT, indent=1))
