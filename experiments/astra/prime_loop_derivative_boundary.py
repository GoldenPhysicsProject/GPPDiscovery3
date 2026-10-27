"""Critical-line prime-power grading: log Euler product versus Weil source.

Tests the important difference between the k-th term of log(zeta) and
the k-th term of -zeta'/zeta at the same prime cutoff.
"""

from __future__ import annotations

import math

import numpy as np


def primes_below(limit: int) -> np.ndarray:
    sieve = np.ones(limit + 1, dtype=np.bool_)
    sieve[:2] = False
    for i in range(2, math.isqrt(limit) + 1):
        if sieve[i]:
            sieve[i * i : limit + 1 : i] = False
    return np.flatnonzero(sieve).astype(np.float64)


def sums(primes: np.ndarray, s: complex, k: int) -> tuple[complex, complex]:
    logp = np.log(primes)
    harmonic = np.exp(-k * s * logp)
    return (complex(np.sum(harmonic) / k), complex(np.sum(logp * harmonic)))


def main() -> None:
    height = 14.134725141734693790457251983562
    cutoffs = [1000, 10000, 100000, 1000000]
    for theta in (0, 0.1):
        print(f"\nshift theta={theta} (lower species s=1/2-theta+i*t)")
        for cutoff in cutoffs:
            primes = primes_below(cutoff)
            s = (0.5 - theta) + 1j * height
            l1, d1 = sums(primes, s, 1)
            l2, d2 = sums(primes, s, 2)
            l3, d3 = sums(primes, s, 3)
            print(
                f"P={cutoff:7d}, log harmonics |k1|={abs(l1):9.5f},"
                f" |k2|={abs(l2):8.5f}, |k3|={abs(l3):8.5f};"
                f" Weil |k2|={abs(d2):8.5f}, |k3|={abs(d3):8.5f}"
            )
            if theta == 0:
                # Integral in the PNT approximation for k=2 Weil boundary.
                # It has a moving, nonconvergent endpoint for every t != 0.
                endpoint = np.exp(-2j * height * math.log(cutoff)) / (-2j * height)
                print(f"  k2 Weil minus oscillatory PNT endpoint: {d2-endpoint:+.6f}")

    print("\nExact coefficient ratios of shifted to original k-th prime mode:")
    for p in (2, 11, 101):
        for k in (1, 2, 3):
            theta = 0.1
            ratio = math.cosh(k * theta * math.log(p))
            print(f"p={p:3}, k={k}: {ratio:.8f}")


if __name__ == "__main__":
    main()
