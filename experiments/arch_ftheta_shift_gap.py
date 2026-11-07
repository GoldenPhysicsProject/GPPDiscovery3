#!/usr/bin/env python3
"""Rigorously discriminate zeta vs F_theta on a zeta-Arch-dark vector.

For the Wave234 *cross-Arch-only* folded kernels, k_theta(w) =
2*cosh(theta*w)*k_zeta(w) with theta=1/4. Thus M_theta=-C_theta
is the sum of Gram squares at lambda_n +/- theta, whereas M_0=-C_0
has only the lambda_n=2n+1/2 squares.

Take coefficients v_j of P_K(z)=prod_{m=0}^{K-2}(z-r_m),
r_m=exp(-(7/10)*lambda_m). Then all first K-1 zeta Laplace
moments vanish. The theta-shifted n=0, minus mode has positive rate
lambda_0-theta=1/4, and evaluates P_K at r_* = exp(-7/40),
strictly larger than every r_m. Its ONE square therefore yields
R_theta(v) >= c_*^2 P_K(r_*)^2 / ||v||^2 > 0.

This script computes a directed 50-decimal-digit rational lower bound
for R_theta(v), and compares it with the *analytic rational upper bound*
on R_zeta(v) from arch_cross_superexp_null.py. Comparisons and signs
are ALL exact (Fraction/integer), not floating-point decisions.

The contrast does not imply an RH proof: full F_theta includes distinct
prime, pole, and diagonal Arch channels with different normalization.
"""

import json
from fractions import Fraction as Q

if __package__:
    try:
        from .arch_hankel_sine2_k6_exact_intervals import (
            SCALE, add, sub, mul, interval, exp_negative, pi_bounds, laplace_sine2)
    except ModuleNotFoundError:
        # Standalone scratch testing, where the established exact source
        # currently lives at the scratch root rather than experiments/.
        from arch_hankel_sine2_k6_exact_intervals import (
            SCALE, add, sub, mul, interval, exp_negative, pi_bounds, laplace_sine2)
    from .arch_cross_superexp_null import exact_upper
else:
    from arch_hankel_sine2_k6_exact_intervals import (
        SCALE, add, sub, mul, interval, exp_negative, pi_bounds, laplace_sine2)
    from arch_cross_superexp_null import exact_upper


def minus(x: tuple[int, int]) -> tuple[int, int]:
    return (-x[1], -x[0])


def theta_mode_rayleigh_lower(k: int) -> Q:
    """Exact rational lower bound using one shifted mode and 50-digit balls."""
    if not 2 <= k <= 15:
        raise ValueError("K must be 2..15 for this directed-ball precision")
    a0, delta, theta = Q(7, 20), Q(7, 10), Q(1, 4)
    lam_star = Q(1, 2) - theta
    star_node = exp_negative(delta * lam_star)

    # Polynomial coefficient vector in ascending monomial order.
    coeff = [interval(1)]
    value = interval(1)
    for m in range(k - 1):
        lam = 2*m + Q(1, 2)
        node = exp_negative(delta * lam)
        diff = sub(star_node, node)
        assert diff[0] > 0, "strict separation of shifted and unshifted rates"
        value = mul(value, diff)
        updated = [(0, 0) for _ in range(len(coeff) + 1)]
        for j, c in enumerate(coeff):
            updated[j] = add(updated[j], minus(mul(node, c)))
            updated[j + 1] = add(updated[j + 1], c)
        coeff = updated

    pi = pi_bounds()
    pi2 = mul(pi, pi)
    cstar = mul(exp_negative(a0 * lam_star), laplace_sine2(lam_star, pi2))
    weighted = mul(cstar, value)
    assert weighted[0] > 0
    numerator = mul(weighted, weighted)
    norm_upper_scaled = sum(mul(c, c)[1] for c in coeff)
    assert norm_upper_scaled > 0 and numerator[0] > 0
    return Q(numerator[0], norm_upper_scaled)


def main() -> None:
    results = []
    for k, minimum_ratio in ((6, 10**14), (10, 10**50)):
        theta_lower = theta_mode_rayleigh_lower(k)
        zeta_upper = exact_upper(k)
        assert theta_lower > minimum_ratio * zeta_upper
        results.append({
            "K": k,
            "zeta_same_vector_analytic_upper_display": f"{float(zeta_upper):.10e}",
            "ftheta_same_vector_interval_lower_display": f"{float(theta_lower):.10e}",
            "certified_ratio_gt": f"1e{len(str(minimum_ratio))-1}",
            "certificate": "R_theta/R_zeta comparison exactly rational, 50-digit outward intervals"
        })
    print(json.dumps({
        "scope": "Wave234 folded cross-Arch only; common zeta-annihilating polynomial vector",
        "theta": "1/4", "results": results,
        "interpretation": "Shifted Laplace frequencies break exact cancellation of the unshifted first K-1 modes.",
        "not_proved": "No full channel positivity and no RH"
    }, indent=2))


if __name__ == "__main__":
    main()
