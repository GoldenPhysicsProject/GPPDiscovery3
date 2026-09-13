"""Truncated prime side of the heat pairing.

Computes
    P_N(t) = 2 * sum_{n<=N} Lambda(n)/sqrt(n) * exp(-(log n)**2 / (4t))
which is the finite piece of GppWeilLadder.primeSide at the heat Gaussian.

This is not the Connes–van Suijlekom matrix Q(c). Do not treat P_N as Q.
"""

from __future__ import annotations

import math
from typing import Iterable


def von_mangoldt(n: int) -> float:
    if n < 2:
        return 0.0
    x = n
    p = 2
    while p * p <= x:
        if x % p == 0:
            while x % p == 0:
                x //= p
            return math.log(p) if x == 1 else 0.0
        p += 1 if p == 2 else 2
    return math.log(n)


def prime_heat(t: float, N: int) -> float:
    if t <= 0:
        raise ValueError("t must be positive")
    acc = 0.0
    for n in range(2, N + 1):
        lam = von_mangoldt(n)
        if lam == 0.0:
            continue
        acc += lam / math.sqrt(n) * math.exp(-(math.log(n) ** 2) / (4.0 * t))
    return 2.0 * acc


def scan(ts: Iterable[float], N: int) -> list[tuple[float, float]]:
    return [(t, prime_heat(t, N)) for t in ts]


if __name__ == "__main__":
    for t, val in scan((0.25, 0.5, 1.0, 2.0), 200):
        print(f"t={t:g}  P_200={val:.12g}")
