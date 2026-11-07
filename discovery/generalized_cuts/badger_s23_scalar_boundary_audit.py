#!/usr/bin/env python3
"""Exact kinematic audit for the s23 bubble in Badger's y,t variables.

Use the four-point bubble choice K1=p2+p3, chi=p2.  In units S1=s23=1,
choose bispinors

    p2 = diag(1,0),  p3 = diag(0,1),  K1 = I.

Then K1^flat=p3 and Badger eq. (39) becomes

    l1 = [[1-y, (y(1-y)-mu2)/t],
          [t,   y                    ]],

with det(l1)=det(l1-K1)=mu2.  For the left tree with scalar momentum -l1,
the remaining scalar propagator is exactly

    D = (-l1+p2)^2-mu2 = -y,

or D=-S1*y after restoring dimensions.

IMPORTANT: this file deliberately makes NO claim that the helicity amplitude itself is
a function of D alone.  A previous version incorrectly inserted a real-frame helicity
representative into Badger's complex cut.  The little-group phase becomes t-dependent
under that continuation.  Using Badger's actual tree formula is required for the bubble
boundary extraction.
"""

from __future__ import annotations

import sympy as sp


y, t, mu2 = sp.symbols("y t mu2", nonzero=True)

p2 = sp.Matrix([[1, 0], [0, 0]])
p3 = sp.Matrix([[0, 0], [0, 1]])
K1 = p2 + p3

l1 = sp.Matrix([
    [1 - y, (y * (1 - y) - mu2) / t],
    [t, y],
])


def sq(p: sp.Matrix):
    return sp.factor(p.det())


def main() -> None:
    assert K1 == sp.eye(2)
    assert sp.simplify(sq(K1) - 1) == 0
    assert sp.simplify(sq(l1) - mu2) == 0
    assert sp.simplify(sq(l1 - K1) - mu2) == 0

    D = sp.factor(sq(-l1 + p2) - mu2)
    assert sp.simplify(D + y) == 0

    print("det(l1) =", sq(l1))
    print("det(l1-K1) =", sq(l1 - K1))
    print("D_left =", D)
    print("PASS: Badger s23 kinematics gives D_left=-y exactly at S1=1")
    print("NOTE: no helicity-amplitude boundary is inferred from D alone")


if __name__ == "__main__":
    main()
