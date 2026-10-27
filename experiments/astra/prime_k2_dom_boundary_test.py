"""Second-prime-harmonic boundary audit on Turbo's actual S1 bump sectors.

The wave-234 bump is sin(pi*u)**2 on [0,1], shifted to
    a_j = .35 + .7*j.
Its correlation G0 has support [-1, 1].  We isolate the k=2 terms of
the EXACT Q (Toeplitz) and C (Hankel) prime kernels, not a surrogate Gram.

Q_ij^prime = -sum_n Lambda(n)/sqrt(n) *
              (G0(log(n)+h*(i-j)) + G0(-log(n)+h*(i-j))).
C_ij^prime = -sum_n Lambda(n)/sqrt(n) *
              G0(log(n)-a_i-a_j-1).

For n=p**2, Lambda(n)/sqrt(n)=(log p)/p.  The finite-support Weil
functional already has the Abel regularized boundary value at sigma=1/2.
Therefore its properly renormalized k=2 matrix is exactly its original
matrix: both P+=Q+C and P-=Q-C have zero matrix correction.

Run: python research/prime_k2_dom_boundary_test.py
"""

from __future__ import annotations

import math


H = 0.7
A0 = 0.35


def g0(x: float) -> float:
    """Autocorrelation of sin^2(pi*u) on [0,1], exact closed expression."""
    y = abs(x)
    if y >= 1.0:
        return 0.0
    return ((1.0 - y) * (2.0 + math.cos(2.0 * math.pi * y)) / 8.0
            + 3.0 * math.sin(2.0 * math.pi * y) / (16.0 * math.pi))


def g0_quadrature(x: float, steps: int = 40000) -> float:
    """Independent midpoint quadrature: only used to verify the formula."""
    shift = abs(x)
    if shift >= 1:
        return 0.0
    width = 1.0 - shift
    return width / steps * sum(
        math.sin(math.pi * (j + .5) * width / steps) ** 2
        * math.sin(math.pi * ((j + .5) * width / steps + shift)) ** 2
        for j in range(steps)
    )


def primes_through(n: int) -> list[int]:
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            sieve[p * p:n + 1:p] = b"\x00" * ((n - p * p) // p + 1)
    return [p for p in range(2, n + 1) if sieve[p]]


def k2_matrices(K: int, prime_cutoff: int, abel_eta: float = 0.0):
    a = [A0 + H * j for j in range(K)]
    q = [[0.0] * K for _ in range(K)]
    c = [[0.0] * K for _ in range(K)]
    for p in primes_through(prime_cutoff):
        y = 2.0 * math.log(p)
        coeff = math.log(p) / p * math.exp(-abel_eta * y)
        for i in range(K):
            for j in range(K):
                q[i][j] -= coeff * (g0(y + H * (i - j))
                                    + g0(-y + H * (i - j)))
                c[i][j] -= coeff * g0(y - a[i] - a[j] - 1.0)
    return q, c


def maxdiff(a, b) -> float:
    return max(abs(x - y) for ra, rb in zip(a, b) for x, y in zip(ra, rb))


def matrix_norm_inf(a) -> float:
    return max(sum(abs(x) for x in row) for row in a)


def k2_hankel_at_center(s: float) -> float:
    """Actual k=2 sewing entry at a_i+a_j+1=s, source sign included."""
    limit = math.ceil(math.exp((s + 1.0) / 2.0))
    return -sum(math.log(p) / p * g0(2.0 * math.log(p) - s)
                for p in primes_through(limit))


def k2_toeplitz_at_lag(v: float) -> float:
    """Actual k=2 diagonal-block coefficient for lag v."""
    limit = math.ceil(math.exp((abs(v) + 1.0) / 2.0))
    return -sum(math.log(p) / p *
                (g0(2.0 * math.log(p) + v)
                 + g0(-2.0 * math.log(p) + v))
                for p in primes_through(limit))


def plusminus(q, c, sign: int):
    return [[x + sign * y for x, y in zip(ra, rb)]
            for ra, rb in zip(q, c)]


def audit() -> None:
    for x in [0, .1, .25, .5, .7, .999]:
        assert abs(g0(x) - g0_quadrature(x)) < 5.e-11, x
    integration_steps = 10000
    correlation_mass = 2.0 / integration_steps * sum(
        g0((j + .5) / integration_steps) for j in range(integration_steps)
    )
    assert abs(correlation_mass - .25) < 1.e-10
    print("G0 exact formula: independent midpoint quadrature PASS")

    print("K  max p^2-support p  norm(Q_k2)  norm(C_k2)  "
          "cutoff delta P+  cutoff delta P-  Abel eta=1e-4 P+ delta")
    for K in range(2, 8):
        amax = A0 + H * (K - 1)
        # Support: cross center a_i+a_j+1 and G0 width 1.
        # Thus log(p^2) < 2*a_max+2; hence p < exp(a_max+1).
        pmax = math.exp(amax + 1.0)
        pfinal = math.ceil(pmax) + 20
        q, c = k2_matrices(K, pfinal)
        qwide, cwide = k2_matrices(K, 2 * pfinal + 111)
        qp, qm = plusminus(q, c, 1), plusminus(q, c, -1)
        qpwide, qmwide = plusminus(qwide, cwide, 1), plusminus(qwide, cwide, -1)
        diff_plus = maxdiff(qp, qpwide)
        diff_minus = maxdiff(qm, qmwide)
        assert diff_plus < 1.e-14 and diff_minus < 1.e-14
        qabel, cabel = k2_matrices(K, pfinal, abel_eta=1.e-4)
        abel_delta = maxdiff(qp, plusminus(qabel, cabel, 1))
        print(f"{K:<2} {pmax:>18.3f} {matrix_norm_inf(q):>11.7f} "
              f"{matrix_norm_inf(c):>11.7f} {diff_plus:>16.2e} "
              f"{diff_minus:>16.2e} {abel_delta:>22.3e}")

    print("\nResult: k=2 has a NONZERO effect on actual S1 parity matrices,")
    print("but its correctly defined Abel/Weil boundary has ZERO extra")
    print("additive correction to either low-K parity pivot matrix.")
    print("A scalar sharp-cutoff oscillation X^(-2it)/(-2it) is not")
    print("a contact term to append to a compactly supported test pairing.")

    # The k=2 boundary DOES have an infinite-volume implication. By PNT,
    # sum_p (log p)/p * G0(2 log p - s) tends (1/2) integral G0 = 1/8.
    # This asymptotic is not a counterterm, but a genuine nondecaying
    # large-separation part of the source, shared by Q and C.
    target = -1.0 / 8.0
    print("\nLarge-shift k=2 source (PNT prediction -1/8):")
    print("center s   C_k2(s)        error from -1/8")
    for center in [6.0, 10.0, 14.0, 18.0, 22.0, 26.0]:
        value = k2_hankel_at_center(center)
        print(f"{center:>7.1f} {value:>14.9f} {value - target:>18.3e}")
        if center == 26.0:
            assert abs(value - target) < .0002
    print("lag v      Q_k2(v)        error from -1/8")
    for lag in [6.0, 10.0, 14.0, 18.0, 22.0]:
        value = k2_toeplitz_at_lag(lag)
        print(f"{lag:>7.1f} {value:>14.9f} {value - target:>18.3e}")


if __name__ == "__main__":
    audit()
