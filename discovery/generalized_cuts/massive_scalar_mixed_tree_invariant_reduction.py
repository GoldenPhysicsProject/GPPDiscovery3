#!/usr/bin/env python3
"""Exact real-frame reduction of the generic mixed-helicity massive-scalar tree.

This audit uses the real centre-of-mass kinematics and the scattering-plane helicity
basis of ``massive_vector_generic_state_sum_symbolic.py``.  In that *fixed real frame*,
with

    s = s23 = 4 E^2,
    mu^2 = E^2 rho^2,
    D  = (q1 + p2)^2,
    D' = (q1 + p3)^2,

one finds exactly

    D + D' = -s,
    A_S(frame) = -2 - 2 D/s - 2 mu^2/D
               =  2 D'/s - 2 mu^2/D.

IMPORTANT: this is not a phase-free Lorentz-scalar representation of a helicity
amplitude.  The chosen real scattering-plane polarization basis fixes a little-group
phase.  Under the complex y,t continuation used in Forde/Badger bubble extraction,
that phase becomes nontrivial and the displayed expression must NOT be substituted for
the actual complex-cut tree.  In particular, a direct evaluation of Badger's mixed
scalar tree on the s23 bubble parametrization gives t-dependent helicity factors.

The identities below remain useful as exact real-kinematics regression tests only.
"""

from __future__ import annotations

import sympy as sp

import massive_vector_generic_state_sum_symbolic as gen
import massive_vector_mhv_state_sum_symbolic as base


def main() -> None:
    ks = gen.generic_kinematics(5)
    g = base.metric(5)
    q1, p2, p3, _q4 = ks

    s = sp.factor(base.mdot(g, p2 + p3, p2 + p3))
    D = sp.factor(base.mdot(g, q1 + p2, q1 + p2))
    Dp = sp.factor(base.mdot(g, q1 + p3, q1 + p3))
    mu2 = sp.factor(base.E**2 * gen.rho**2)

    AS_pm = sp.factor(gen.scalar_tree(+1, -1))
    AS_mp = sp.factor(gen.scalar_tree(-1, +1))

    assert sp.simplify(s - 4 * base.E**2) == 0
    assert sp.simplify(D + Dp + s) == 0
    assert sp.simplify(AS_pm - AS_mp) == 0

    frame_form = sp.factor(-2 - 2 * D / s - 2 * mu2 / D)
    complementary_form = sp.factor(2 * Dp / s - 2 * mu2 / D)

    assert sp.simplify(AS_pm - frame_form) == 0
    assert sp.simplify(frame_form - complementary_form) == 0

    print("s23 =", s)
    print("D =", D)
    print("D' =", Dp)
    print("D + D' =", sp.factor(D + Dp))
    print("mu^2 =", mu2)
    print("A_S in fixed real helicity frame =", AS_pm)
    print("frame reduction =", frame_form)
    print("PASS: exact real-frame identity; NOT valid as a phase-free complex-cut formula")


if __name__ == "__main__":
    main()
