"""Exact pole/first-prime PNT mean-field Ward cancellation on Turbo's S1 sector.

This constructs an actual OS-reflection-positive mean-field covariance
K(u,v)=exp(-abs(u-v)/2) on the doubled prime orientation line.
It does NOT include the signed Chebyshev error measure or k>=2 harmonics,
and therefore is NOT a RH proof.  The identity is exact algebra.

Run: python research/prime_pnt_os_ward_background.py
"""

import math

import numpy as np

from prime_k2_dom_boundary_test import (
    A0, H, g0, bump_laplace, correlation_laplace,
)


def midpoint(f, lo, hi, n=12000):
    width = (hi - lo) / n
    return width * sum(f(lo + (j + .5) * width) for j in range(n))


def pole_q(v):
    m = correlation_laplace(.5)
    return m * (math.exp(v / 2) + math.exp(-v / 2))


def pnt_first_prime_q(v):
    """The PNT main prime-side source, with the original minus sign."""
    return -midpoint(
        lambda u: math.exp(u/2) * (g0(u+v) + g0(v-u)),
        0.0, abs(v) + 1.0,
    )


def massive_covariance_smeared(v):
    return midpoint(lambda w: g0(w) * math.exp(-abs(v-w)/2), -1., 1.)


def pole_c(a, b):
    plus = bump_laplace(.5)
    minus = bump_laplace(-.5)
    return (math.exp((a+b)/2) * plus**2
            + math.exp(-(a+b)/2) * minus**2)


def pnt_first_prime_c(a, b):
    s = a + b + 1.0
    return -math.exp(s/2) * correlation_laplace(.5)


def reflected_covariance(a, b):
    return (math.exp(-(a+b)/2) * bump_laplace(-.5)**2)


def run():
    # Independent quadrature of pole minus unsmeared leading prime continuum.
    print("lag v    pole+PNT prime Q      massive Gram Q        delta")
    for v in [0., .35, .7, 1.4, 2.1, 3.5, 4.2]:
        direct = pole_q(v) + pnt_first_prime_q(v)
        gram = massive_covariance_smeared(v)
        difference = direct - gram
        print(f"{v:>5.2f} {direct:>20.12f} {gram:>20.12f} {difference:>12.2e}")
        assert abs(difference) < 2.e-8

    print("\ncross a,b    pole+PNT prime C     reflected OS C       delta")
    for a,b in [(A0,A0), (A0,A0+H), (A0+2*H,A0+5*H)]:
        direct = pole_c(a,b) + pnt_first_prime_c(a,b)
        os_value = reflected_covariance(a,b)
        difference = direct - os_value
        print(f"{a:>4.2f},{b:<4.2f} {direct:>20.12f} {os_value:>20.12f} {difference:>12.2e}")
        assert abs(difference) < 2.e-14

    # Positive semidefinite both parities because this is the restriction
    # of a positive covariance exp(-|u-v|/2) on the doubled full line.
    print("\nK   min eig(Q+C)     min eig(Q-C)     min eig(Q)     OS cross rank")
    for k in range(2, 8):
        a = [A0 + H*j for j in range(k)]
        q = np.array([[massive_covariance_smeared(a[i]-a[j])
                       for j in range(k)] for i in range(k)])
        c = np.array([[reflected_covariance(a[i],a[j])
                       for j in range(k)] for i in range(k)])
        lp = np.linalg.eigvalsh(q+c)[0]
        lm = np.linalg.eigvalsh(q-c)[0]
        lq = np.linalg.eigvalsh(q)[0]
        cross_rank = np.linalg.matrix_rank(c, tol=1.e-11)
        print(f"{k}  {lp:>14.8e} {lm:>14.8e} {lq:>14.8e} {cross_rank:>10d}")
        assert lp > -1.e-10 and lm > -1.e-10 and cross_rank == 1

    print("\nExact underlying algebra: Pole+PNT(mu1^bulk) =")
    print("  Q_ij = integral G0(w) exp(-|a_i-a_j-w|/2) dw")
    print("  C_ij = exp(-(a_i+a_j)/2) S(-1/2)^2, a PSD OS rank-one cross.")
    print("ALL tests passed. The actual first-prime Chebyshev error remains open.")


if __name__ == "__main__":
    run()
