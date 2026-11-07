#!/usr/bin/env python3
"""36 rigorous intervals for the Wave234 sine-squared CROSS-arch matrix.

Only stdlib integer/Fraction arithmetic is used: fixed-point intervals
are rounded outward, pi is enclosed by Machin's alternating arctangent
series, and exp(-x) by positive Taylor bounds after dyadic scaling.

This is the CROSS-ARCH block of Wave234's S1 construction, not necessarily
the distinct rounded Phat6 matrix. It does not enclose the diagonal arch,
prime/pole channels, or an all-K RH certificate.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import factorial

SCALE = 10**50
TAIL_MASS = Fraction(800, 32**6 * 5_000_000_000)


def ceildiv(a: int, b: int) -> int:
    assert b > 0
    return -((-a) // b)


def interval(q: Fraction | int) -> tuple[int, int]:
    q = Fraction(q)
    return (q.numerator * SCALE // q.denominator,
            ceildiv(q.numerator * SCALE, q.denominator))


def add(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return (x[0] + y[0], x[1] + y[1])


def sub(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return (x[0] - y[1], x[1] - y[0])


def mul(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    pp = (x[0]*y[0], x[0]*y[1], x[1]*y[0], x[1]*y[1])
    return (min(pp) // SCALE, ceildiv(max(pp), SCALE))


def div(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    assert 0 < y[0], "positive denominator interval required"
    pp = (Fraction(x[0]*SCALE, y[0]), Fraction(x[0]*SCALE, y[1]),
          Fraction(x[1]*SCALE, y[0]), Fraction(x[1]*SCALE, y[1]))
    return (int(min(pp).__floor__()), int(max(pp).__ceil__()))


def atan_sum(x: Fraction, degree: int) -> Fraction:
    return sum(((-1)**k * x**(2*k+1) / (2*k+1)
                for k in range(degree+1)), Fraction(0))


def pi_bounds() -> tuple[int, int]:
    # Machin's identity pi=16atan(1/5)-4atan(1/239).
    # Alternating-series lower bounds at ODD degrees, upper at EVEN.
    a5_low, a5_hi = atan_sum(Fraction(1, 5), 41), atan_sum(Fraction(1, 5), 40)
    a239_low, a239_hi = atan_sum(Fraction(1, 239), 11), atan_sum(Fraction(1, 239), 10)
    return (interval(16*a5_low-4*a239_hi)[0],
            interval(16*a5_hi-4*a239_low)[1])


def exp_negative(x: Fraction) -> tuple[int, int]:
    """Mathematically enclosed e^(-x), x>=0, by rational estimates."""
    assert x >= 0
    if not x:
        return (SCALE, SCALE)
    halves = 0
    while x > 1:
        x /= 2
        halves += 1
    n = 70
    s = sum((x**k/factorial(k) for k in range(n+1)), Fraction(0))
    next_term = x**(n+1)/factorial(n+1)
    upper_exp = s + next_term/(1 - x/(n+2))
    out = (interval(1/upper_exp)[0], interval(1/s)[1])
    for _ in range(halves):
        out = mul(out, out)
    return out


def laplace_sine2(lam: Fraction, pi2: tuple[int, int]) -> tuple[int, int]:
    # S(-lam)=2*pi²*(1-e^-lam)/(lam*(lam²+4*pi²)); all factors >0.
    one_minus = sub(interval(1), exp_negative(lam))
    numerator = mul(interval(2), pi2)
    denominator = mul(interval(lam), add(interval(lam**2), mul(interval(4), pi2)))
    return mul(div(numerator, denominator), one_minus)


def cert_tail() -> None:
    # n>=16: lam>=65/2>32 and min (a_j+a_k)=7/10.
    # exp((65/2)*(7/10))=exp(91/4) exceeds 5,000,000,000,
    # using only the 26th positive rational Taylor partial sum.
    x = Fraction(91, 4)
    t = sum((x**k/factorial(k) for k in range(27)), Fraction(0))
    assert t > 5_000_000_000
    assert 1 + Fraction(7, 5) + Fraction(7, 5)**2/2 > 2
    assert TAIL_MASS < Fraction(2, 10**16)


def main() -> None:
    cert_tail()
    pi = pi_bounds()
    pi2 = mul(pi, pi)
    lam = [Fraction(4*n+1, 2) for n in range(16)]
    pos = [Fraction(7, 20) + Fraction(7*j, 10) for j in range(6)]
    s = [laplace_sine2(t, pi2) for t in lam]
    ell = [[mul(exp_negative(t * a), sn) for a in pos]
           for t, sn in zip(lam, s)]
    tail_hi = interval(TAIL_MASS)[1]
    matrix = []
    worst_width = 0
    for j in range(6):
        row = []
        for k in range(6):
            head = (0, 0)
            for n in range(16):
                head = add(head, mul(ell[n][j], ell[n][k]))
            lo, hi = (-head[1] - tail_hi, -head[0])
            assert lo <= hi < 0
            worst_width = max(worst_width, hi - lo)
            row.append({"lower_scaled":str(lo), "upper_scaled":str(hi)})
        matrix.append(row)
    assert all(matrix[j][k] == matrix[k][j] for j in range(6) for k in range(6))
    # Comparison to 2e-16 is exact. Each scaled width is <SCALE/5e15.
    assert worst_width < SCALE // (5*10**15)
    print(json.dumps({"scale":str(SCALE),
                      "description":"Wave234 sine2 6x6 CROSS-arch only, exact outward integer intervals",
                      "terms_used":16,
                      "uniform_tail_bound":"800/(32^6*5000000000) < 2e-16",
                      "worst_interval_width_scaled":str(worst_width),
                      "matrix":matrix}, indent=2))


if __name__ == "__main__":
    main()
