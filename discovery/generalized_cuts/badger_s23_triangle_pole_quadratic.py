#!/usr/bin/env python3
"""Exact quadratic structure of the s23 triangle-subtraction pole.

Companion to ``badger_s23_mhv_pure_bubble_boundary.py``.  In the same explicit
S1=s23=1, Q=1 bispinor frame, the second uncut scalar propagator appearing in
the right mixed-helicity tree has denominator proportional to

    P(y) = u y^2 + (t(1-u^2)-u)y + u mu2 - u t^2 + u^2 t.

The two triangle-subtraction solutions are the two roots of P.  Their
discriminant collapses exactly to

    Delta_y = (t(1+u^2)-u)^2 - 4 mu2 u^2.

This is the square-root structure entering Badger's y_+/- subtraction.  The
Vieta relations provide a branch-free way to sum the two solutions, and their
large-t branches begin as y_+ ~ u t and y_- ~ -t/u.

The file identifies the subtraction geometry only; it does not yet claim the
normalized triangle-subtraction coefficient.
"""

from __future__ import annotations

import sympy as sp

u, y, t, mu2 = sp.symbols("u y t mu2", nonzero=True)

P = sp.expand(
    u * y**2 + (t * (1 - u**2) - u) * y
    + u * mu2 - u * t**2 + u**2 * t
)


def main() -> None:
    disc = sp.factor(sp.discriminant(P, y))
    target_disc = sp.expand((t * (1 + u**2) - u)**2 - 4 * mu2 * u**2)
    assert sp.simplify(disc - target_disc) == 0

    B = t * (1 - u**2) - u
    C = u * mu2 - u * t**2 + u**2 * t
    root_sum = sp.factor(-B / u)
    root_prod = sp.factor(C / u)

    assert sp.simplify(root_sum - (t * (u**2 - 1) + u) / u) == 0
    assert sp.simplify(root_prod - (mu2 - t**2 + u * t)) == 0

    # Branch-leading coefficients obtained directly from the exact quadratic:
    # writing y ~ a t, the t^2 coefficient is u a^2 + (1-u^2)a - u.
    a = sp.symbols("a")
    leading = sp.factor(u * a**2 + (1 - u**2) * a - u)
    assert sp.factor(leading - (a - u) * (u * a + 1)) == 0

    print("P(y) =", sp.factor(P))
    print("Delta_y =", disc)
    print("y_+ + y_- =", root_sum)
    print("y_+ y_- =", root_prod)
    print("leading branch equation =", leading)
    print("large-t slopes: a_+=u, a_-=-1/u")
    print("PASS: triangle subtraction roots have exact branch-free quadratic structure")


if __name__ == "__main__":
    main()
