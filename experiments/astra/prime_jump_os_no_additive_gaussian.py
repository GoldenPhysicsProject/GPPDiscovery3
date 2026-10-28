"""Exact high-frequency bound against naive additive gamma + prime jumps.

The proof is analytic and given in the companion note. This executable
prints certified upper bounds for kinetic derivatives using g'(q)<=6/q,
and numerical examples. No zeta zeros, conjectures, or external packages.
"""

from math import log, pi, sin, sqrt


def gamma_derivative(q, terms=100_000):
    """Finite positive partial sum of exact gamma complete-Bernstein derivative."""
    return sum(2 * (2 * n + 0.5) / (((2 * n + 0.5) ** 2 + q) ** 2)
               for n in range(terms))


def gamma_derivative_upper(q):
    """Rigorous simple uniform bound for q>=1, proved by splitting a_n<=sqrt(q)."""
    assert q >= 1
    return 6.0 / q


def jump_derivative(t, primes):
    """d/dq of sum_p [(log p)/sqrt(p)]*(1-cos(t log p)), q=t*t."""
    return sum((log(p) / sqrt(p)) * log(p) * sin(t * log(p)) / (2 * t)
               for p in primes)


def main():
    c = log(2) / sqrt(2)
    ell = log(2)
    examples = []
    # sin(ell*t_n)=-1 exactly by construction (up to floating roundoff).
    for n in (5, 10, 20, 40):
        t = (4 * n + 3) * pi / (2 * ell)
        bound = gamma_derivative_upper(t * t) - c * ell / (2 * t)
        assert bound < 0, (n, t, bound)
        actual = gamma_derivative(t * t, terms=10_000) - c * ell / (2 * t)
        assert actual < bound
        examples.append((n, round(t, 3), round(bound, 8)))

    # Pair frequencies are incommensurable: no phase synchronisation assumed.
    # A negative bound from the exact gamma inequality suffices for certification.
    pair_certificates = []
    for tenth in range(500, 10001):
        t = tenth / 10
        signed_upper = gamma_derivative_upper(t * t) + jump_derivative(t, (2, 3))
        if signed_upper < 0:
            pair_certificates.append((round(t, 2), round(signed_upper, 8)))
            if len(pair_certificates) >= 5:
                break
    assert len(pair_certificates) == 5

    print("Analytic uniform gamma derivative bound: g'(q) <= 6/q (q>=1)")
    print("Single prime p=2 certified f'(t_n^2)<0 (n,t,upper):", examples)
    print("Two primes p=2,3 certified negative (t,upper):", pair_certificates)
    print("PASS: additive prime jumps cannot preserve gamma Gaussian OS for all scales")


if __name__ == "__main__":
    main()
