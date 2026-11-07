#!/usr/bin/env python3
"""Exact convention-free algebra behind the closed s23 scalar bubble.

This file deliberately strips away all cut-normalization and scalar-multiplicity
conventions.  It certifies only the polynomial/rational identities that may be
safely promoted to a proof assistant.

In the rational frame used by the companion Badger audits,

    s23 = 1,
    s12 = -u^2/(1+u^2),
    Q   = 1.

The one-flow triangle-subtraction coefficient found from the explicit triple cut is

    c_one = -i (5 u^2 + 3)/(3(1+u^2)).

If an external convention supplies multiplicity two (as for one complex scalar),
then 2*c_one equals exactly the phase-normal invariant target

    (2 i / 3) (2 s12 - 3 s23)

in this frame.

The same file also certifies the triangle-root discriminant identity used to remove
the y_± square root by symmetric summation.
"""

from __future__ import annotations

import sympy as sp

u, t, mu2 = sp.symbols("u t mu2", real=True)
I = sp.I


def main() -> None:
    # Triangle-pole quadratic P(y)=u y^2 + B y + C.
    B = t * (1 - u**2) - u
    C = u * mu2 - u * t**2 + u**2 * t
    discriminant = sp.expand(B**2 - 4 * u * C)
    discriminant_target = sp.expand((t * (1 + u**2) - u)**2 - 4 * mu2 * u**2)
    assert sp.expand(discriminant - discriminant_target) == 0

    # Frame invariant and one-flow result.
    s23 = sp.Integer(1)
    s12 = -u**2 / (1 + u**2)
    c_one = -I * (5 * u**2 + 3) / (3 * (1 + u**2))
    target_frame = sp.Rational(2, 3) * I * (2 * s12 - 3 * s23)

    # This equality is purely algebraic.  The factor 2 is NOT derived here; it is
    # an explicit multiplicity supplied by the complex-scalar convention.
    assert sp.factor(2 * c_one - target_frame) == 0

    # Real coefficient version, ideal for Lean formalization without complex
    # arithmetic obscuring the rational identity.
    real_one = -(5 * u**2 + 3) / (3 * (1 + u**2))
    real_target = sp.Rational(2, 3) * (2 * s12 - 3 * s23)
    assert sp.factor(2 * real_one - real_target) == 0

    # 1+u^2 is strictly positive on the real frame, so no hidden exceptional
    # kinematic point was introduced by the rational restoration.
    # SymPy cannot prove positivity for a symbolic real u by simplification alone;
    # encode the exact sum-of-squares decomposition instead.
    denom_minus_one = sp.expand((1 + u**2) - 1)
    assert denom_minus_one == u**2

    print("Delta_y =", sp.factor(discriminant))
    print("s12 =", s12)
    print("one-flow coefficient =", sp.factor(c_one))
    print("twice one-flow =", sp.factor(2 * c_one))
    print("frame target =", sp.factor(target_frame))
    print("PASS: discriminant and frame restoration identities are exact")
    print("NOTE: multiplicity two remains an explicit complex-scalar convention")


if __name__ == "__main__":
    main()
