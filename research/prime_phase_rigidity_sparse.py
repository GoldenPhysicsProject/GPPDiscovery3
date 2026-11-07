"""Independent-prime clock diagnostics, with no RH claims.

The exact mathematical proofs are in PRIME_PHASE_RIGIDITY_SPARSE.md.
This is a standard-library numerical witness and finite-factor checker.
"""

from __future__ import annotations

import cmath
from decimal import Decimal, localcontext
from fractions import Fraction


def convergents(x: Decimal, count: int) -> list[Fraction]:
    """Continued-fraction convergents of a positive high-precision Decimal."""
    if x <= 0 or count < 0:
        raise ValueError("positive x and nonnegative count required")
    out: list[Fraction] = []
    a0, a1 = 0, 1
    b0, b1 = 1, 0
    y = x
    for _ in range(count):
        a = int(y)
        num, den = a * a1 + a0, a * b1 + b0
        out.append(Fraction(num, den))
        a0, a1, b0, b1 = a1, num, b1, den
        rem = y - a
        if rem == 0:
            break
        y = 1 / rem
    return out


def real_dilation_discontinuity(precision: int = 90, count: int = 11) -> list[dict]:
    """Witness p^a/q^b -> 1, while frozen-p clock diverges.

    Here alpha_2=0 and alpha_3=1.  A convergent b/a to log(2)/log(3)
    makes a*log(2)-b*log(3) small, but the altered clock is -b*log(3).
    Returned decimal approximations illustrate a separate exact proof.
    """
    if precision < 30:
        raise ValueError("precision too low for the demonstration")
    with localcontext() as ctx:
        ctx.prec = precision
        log2, log3 = Decimal(2).ln(), Decimal(3).ln()
        rows = []
        for f in convergents(log2 / log3, count):
            b, a = f.numerator, f.denominator
            original = Decimal(a) * log2 - Decimal(b) * log3
            slowed = -Decimal(b) * log3
            rows.append({
                "a": a, "b": b,
                "log_n": str(+original),
                "altered_log_n": str(+slowed),
            })
        return rows


def euler_correction(s: complex, speeds: dict[int, float]) -> complex:
    """Finite M(s)=prod_p (1-p^-s)/(1-p^-s'_p)."""
    if s.real <= 0:
        raise ValueError("Re(s) must be positive")
    result = 1 + 0j
    for p, alpha in speeds.items():
        if p < 2 or not (0 <= alpha <= 1):
            raise ValueError("primes >= 2 and rates in [0,1] required")
        # Caller supplies actual primes; this is not a primality verifier.
        sp = 0.5 + alpha * (s - 0.5)
        result *= (1 - cmath.exp(-s * cmath.log(p))) / (
            1 - cmath.exp(-sp * cmath.log(p))
        )
    return result


def reflection_completion(s: complex, speeds: dict[int, float]) -> complex:
    """Finite H(s)=M(s)M(1-s), holomorphic on the critical strip."""
    if not (0 < s.real < 1):
        raise ValueError("critical strip required")
    return euler_correction(s, speeds) * euler_correction(1-s, speeds)


def upper_bound_on_line(speeds: dict[int, float]) -> float:
    """Uniform-in-t bound H(1/2+it) <= this (finite product)."""
    bound = 1.0
    for p in speeds:
        r = p ** (-0.5)
        bound *= ((1 + r) / (1 - r)) ** 2
    return bound


if __name__ == "__main__":
    print("Continuous real-dilation clock obstruction (2 frozen, 3 unchanged)")
    for row in real_dilation_discontinuity():
        print(f"2^{row['a']}/3^{row['b']}  "
              f"log(n)={Decimal(row['log_n']):+.5E}  "
              f"altered={Decimal(row['altered_log_n']):+.5E}")
    speeds = {2: 0, 3: 0.25, 5: 0.8, 7: 1}
    print("FE-completed finite factor at critical height 14.134725...")
    z = reflection_completion(complex(0.5, 14.134725141734693), speeds)
    print(f"H={z.real:.12g}{z.imag:+.2g}j, uniform upper={upper_bound_on_line(speeds):.6g}")
