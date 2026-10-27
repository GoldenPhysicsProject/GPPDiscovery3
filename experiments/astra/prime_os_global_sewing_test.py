"""Numerical audits of doubled-prime OS reflection and the scaling window.

Exact local algebra is recorded in the accompanying dated research note.
This program tests the matrices, not a global RH claim.
"""

from __future__ import annotations

import math
import cmath

import numpy as np
import mpmath as mp


def primes_below(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for n in range(2, math.isqrt(limit) + 1):
        if sieve[n]:
            sieve[n * n : limit + 1 : n] = b"\x00" * (
                (limit - n * n) // n + 1
            )
    return [i for i, is_prime in enumerate(sieve) if is_prime]


def doubled_kernels(p: int, sigmas: list[float]) -> tuple[np.ndarray, np.ndarray]:
    """Return raw swap reflection and ordinary positive doubled-field Gram."""
    r = 1 / p
    x = np.array([(sigma - 0.5) * math.log(p) for sigma in sigmas])
    d = np.sqrt((1 - r * np.exp(2 * x)) * (1 - r * np.exp(-2 * x)))
    difference = x[:, None] - x[None, :]
    total = x[:, None] + x[None, :]
    os_raw = np.outer(d, d) / (
        (1 - r * np.exp(difference)) * (1 - r * np.exp(-difference))
    )
    gram_positive = np.outer(d, d) / (
        (1 - r * np.exp(total)) * (1 - r * np.exp(-total))
    )
    return os_raw, gram_positive


def crossover_log_overlap(primes: list[int], c: float, d: float) -> complex:
    """Sigma/t complex shift z=(c+id)/log P, relative to the midpoint."""
    L = math.log(primes[-1])
    value = 0j
    for p in primes:
        u = math.log(p) / L
        r = 1 / p
        q0 = math.sqrt(r)
        q1 = q0 * cmath.exp(-(c + 1j * d) * u)
        value += 0.5 * math.log(1 - abs(q0) ** 2)
        value += 0.5 * math.log(1 - abs(q1) ** 2)
        value -= cmath.log(1 - np.conj(q0) * q1)
    return value


def crossover_integral(c: float, d: float) -> complex:
    """PNT limiting log coherent inner product, including phase."""
    mp.mp.dps = 35
    f = lambda u: (
        mp.exp(-(c + 1j * d) * u)
        - mp.mpf("0.5")
        - mp.mpf("0.5") * mp.exp(-2 * c * u)
    ) / u
    # The removable endpoint value is -(1j*d), not numerical 0/0.
    return complex(mp.quad(f, [0, 1]))


def mirror_log_overlap(primes: list[int], c: float) -> float:
    L = math.log(primes[-1])
    answer = 0.0
    for p in primes:
        sigma = 0.5 + c / L
        r = 1 / p
        a = p ** (-sigma)
        b = p ** (-(1 - sigma))
        answer += (math.log1p(-a * a) + math.log1p(-b * b)) / 2
        answer -= math.log1p(-r)
    return answer


def mirror_limit(c: float) -> float:
    mp.mp.dps = 35
    return float(-mp.quad(lambda u: (mp.cosh(2 * c * u) - 1) / u, [0, 1]))


def main() -> None:
    sigmas = [0.52, 0.63, 0.81]
    for p in (2, 3, 5, 17):
        raw, corrected = doubled_kernels(p, sigmas)
        raw_eig = np.linalg.eigvalsh(raw)
        corrected_eig = np.linalg.eigvalsh(corrected)
        print(f"p={p:2}: raw J-swap eigs {raw_eig}, ordinary Gram eigs {corrected_eig}")
        assert raw_eig[0] < -1e-13
        assert corrected_eig[0] > -1e-12
        # A 2x2 minor is already negative, without a numerical eigen solver.
        assert np.linalg.det(raw[:2, :2]) < 0
        assert abs(corrected[0, 0] - 1.0) < 1e-12

    for limit in (1000, 10000, 100000, 1000000):
        primes = primes_below(limit)
        c, d = 0.8, 0.35
        overlap = crossover_log_overlap(primes, c, d)
        mir = mirror_log_overlap(primes, c)
        print(f"P={limit:7}: crossover logK {overlap.real:+.8f}{overlap.imag:+.8f}i, mirror logF {mir:+.8f}")
    print("Predicted crossover logK:", crossover_integral(c, d))
    print("Predicted mirror logF:", mirror_limit(c))

    # Exact structural identity from the doubly oriented formula.
    for p in (2, 3, 5):
        r = 1 / p
        for x in (0.04, 0.12, -0.12):
            h = 1 / ((1 - r * math.exp(x)) * (1 - r * math.exp(-x)))
            s = sum(r ** abs(k) * math.exp(k * x) for k in range(-120, 121)) / (1 - r * r)
            assert abs(h - s) < 1e-10
    print("Bilateral hyperbolic mode expansion: PASS")


if __name__ == "__main__":
    main()
