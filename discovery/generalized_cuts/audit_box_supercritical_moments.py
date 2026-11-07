#!/usr/bin/env python3
"""Exact audit of the all-order supercritical massless-box moment formula.

No floating-point arithmetic is used.  For each m we compare

  sum_j binomial(m,j) s^j t^(m-j)
        * (j!)^2 ((m-j)!)^2 / (2m+3)!

against

  m!/(2m+3)! * sum_j j!(m-j)! s^j t^(m-j),

then multiply by the n=4 supercritical coefficient
(-1)^(m+1)(m+1).
"""

from sympy import symbols, factorial, binomial, expand, simplify

s, t = symbols("s t")


def expanded_dirichlet_moment(m: int):
    denom = factorial(2 * m + 3)
    return expand(
        sum(
            binomial(m, j)
            * s**j
            * t ** (m - j)
            * factorial(j) ** 2
            * factorial(m - j) ** 2
            / denom
            for j in range(m + 1)
        )
    )


def closed_dirichlet_moment(m: int):
    return expand(
        factorial(m)
        / factorial(2 * m + 3)
        * sum(
            factorial(j) * factorial(m - j) * s**j * t ** (m - j)
            for j in range(m + 1)
        )
    )


def box_residue(m: int):
    return expand((-1) ** (m + 1) * (m + 1) * closed_dirichlet_moment(m))


def main():
    for m in range(9):
        lhs = expanded_dirichlet_moment(m)
        rhs = closed_dirichlet_moment(m)
        assert simplify(lhs - rhs) == 0
        assert simplify(rhs - rhs.xreplace({s: t, t: s})) == 0
        print(f"m={m}: moment = {rhs}")
        print(f"     residue = {box_residue(m)}")

    assert simplify(box_residue(0) + 1 / 6) == 0
    assert simplify(box_residue(1) - (s + t) / 60) == 0
    assert simplify(
        box_residue(2) + (2 * s**2 + s * t + 2 * t**2) / 840
    ) == 0
    assert simplify(
        box_residue(3)
        - (s**3 / 2520 + (s**2 * t + s * t**2) / 7560 + t**3 / 2520)
    ) == 0

    print("All exact audits passed.")


if __name__ == "__main__":
    main()
