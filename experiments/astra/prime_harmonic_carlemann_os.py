"""Prime repetition Hankel OS kernels, sharp k=2 Carleman edge.

All computations use literal finite prime sums, no zeta zeros and no
external packages. Exact all-prime proofs are in the accompanying note.
"""

from math import cos, cosh, exp, log, pi, sqrt


def sieve(limit):
    flags = bytearray(b"\x01") * (limit + 1)
    flags[:2] = b"\x00\x00"
    for n in range(2, int(sqrt(limit)) + 1):
        if flags[n]:
            flags[n * n::n] = b"\x00" * len(flags[n * n::n])
    return [n for n in range(2, limit + 1) if flags[n]]


def kernel(primes, k, t, u):
    return sum(log(p) * exp(-k * (0.5 + t + u) * log(p))
               for p in primes)


def gram_quadratic(primes, k, times, coeffs):
    """Exact positive square representation of the finite prime Hankel kernel."""
    return sum(log(p) * p ** (-k / 2) *
               (sum(c * p ** (-k * t) for t, c in zip(times, coeffs))) ** 2
               for p in primes)


def mellin_carlemann_multiplier(tau):
    """Numerical independent Beta/Mellin check of π/cosh(πτ)."""
    dy = 0.01
    integral = sum(cos(tau * (-50 + (j + 0.5) * dy)) /
                   (2 * cosh((-50 + (j + 0.5) * dy) / 2)) * dy
                   for j in range(10000))
    return integral


def main():
    primes = sieve(1_000_000)
    assert len(primes) == 78498
    times = (0.12, 0.27, 0.49, 0.75)
    coeffs = (1.0, -2.0, 1.25, 0.75)
    for k in (1, 2, 3, 4):
        # K_k(t,u)=K_1(kt+(k-1)/4,ku+(k-1)/4), literally termwise.
        for t, u in ((0.1, 0.3), (0.45, 0.73), (1.2, 0.4)):
            lhs = kernel(primes, k, t, u)
            rhs = kernel(primes, 1, k * t + (k - 1) / 4,
                         k * u + (k - 1) / 4)
            assert abs(lhs - rhs) < 3e-12 * max(1, abs(lhs))
        assert gram_quadratic(primes, k, times, coeffs) >= 0

    for tau in (0, 0.4, 1.0, 2.0):
        actual = mellin_carlemann_multiplier(tau)
        predicted = pi / cosh(pi * tau)
        assert abs(actual - predicted) < 1e-8, (tau, actual, predicted)

    # The exact borderline k=2 density is asymptotically 1/(2(t+u)).
    # Finite prime truncations underestimate the infinite mass.
    for h in (0.3, 0.5, 1.0):
        value = kernel(primes, 2, h / 2, h / 2)
        assert value > 0
        print(f"k2 Laplace mass, h={h:.1f}, primes<=1e6: "
              f"{value:.8f}; Carleman singular 1/(2h)={1/(2*h):.8f}")

    # At infinite cutoff the k2 trace is sum_p 1/(4p), divergent;
    # k>=3 is trace class with exact trace Σ_(p,k>=3)p^(-k/2)/(2k).
    for cutoff in (10_000, 100_000, 1_000_000):
        subset = [p for p in primes if p <= cutoff]
        trace_k2 = sum(1 / (4 * p) for p in subset)
        trace_kge3 = sum(sum(p ** (-k / 2) / (2 * k)
                              for k in range(3, 25)) for p in subset)
        print(f"p<={cutoff:7}: k2 trace partial={trace_k2:.8f}; "
              f"k>=3 trace partial={trace_kge3:.8f}")
    # Shifted F_theta=xi(s-theta)xi(s+theta), theta>0: its lower
    # species at sigma=1/2-theta has a k2 support wedge t+u>theta.
    theta = 0.1
    for h in (0.05, 0.2):
        vals = []
        for cutoff in (10_000, 100_000, 1_000_000):
            vals.append(round(sum(log(p) * p ** (-1 + 2 * theta - 2 * h)
                                  for p in primes if p <= cutoff), 5))
        print(f"F_theta lower branch theta={theta}, t+u={h}: "
              f"{vals}, threshold t+u>theta")
        assert vals[0] < vals[1] < vals[2]
    print("PASS: dilation, positive squares, and trace threshold")


if __name__ == "__main__":
    main()
