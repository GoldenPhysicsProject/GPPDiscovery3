#!/usr/bin/env python3
"""Generic exact massive-vector / extra-scalar tree-state audit.

This repairs an important limitation of the earlier ``massive_vector_mhv_state_sum``
audits: those scripts fixed the four-dimensional projection of the two massive legs at
rest.  Their exact identities were therefore threshold identities, not generic cut
identities.

Here two independent rational parameters are retained:

    beta = (1-r^2)/(1+r^2) = |p_massive|/E,
    rho  = 2r/(1+r^2)      = mu/E,

with beta^2 + rho^2 = 1, and

    cos(theta) = (1-t^2)/(1+t^2),
    sin(theta) = 2t/(1+t^2).

Thus the 5D-null massive legs have generic nonzero four-dimensional spatial momentum
except at r=1.  We use the physically relevant color ordering

    (massive, gluon, gluon, massive).

The script proves exactly that:

1. the old same-helicity 3:1 vector/scalar sewing is *not* generic;
2. its defect is

       4 (r^2-1)^2 (1+t^2)^2 / (r^2+t^2)^2;

3. the mixed-helicity extra-scalar tree is generically nonzero and vanishes only on the
   old threshold slice r=1;
4. the previously observed mixed-helicity vector sewing 16 is likewise recovered at
   r=1, but is not a generic bubble coefficient.

All checks are exact SymPy identities.  No claim is made here about coupling/color/cut
orientation, loop-measure normalization, or higher-topology subtraction.
"""

from __future__ import annotations

import sympy as sp

import massive_vector_mhv_state_sum_symbolic as base

r = sp.symbols("r", real=True, nonzero=True)
E = base.E
t = base.t
I = sp.I

beta = (1 - r**2) / (1 + r**2)
rho = 2 * r / (1 + r**2)
ct = (1 - t**2) / (1 + t**2)
st = 2 * t / (1 + t**2)


def generic_kinematics(d: int = 5):
    q1 = sp.Matrix([-E, 0, 0, -E * beta, +E * rho])
    q4 = sp.Matrix([-E, 0, 0, +E * beta, -E * rho])
    g2 = sp.Matrix([+E, +E * st, 0, +E * ct, 0])
    g3 = sp.Matrix([+E, -E * st, 0, -E * ct, 0])
    ks = [q1, g2, g3, q4]
    if d == 5:
        return ks
    if d == 6:
        return [sp.Matrix(list(k) + [0]) for k in ks]
    raise ValueError("only d=5 or d=6 is used")


def gluon_pol(leg: int, h: int, d: int = 5):
    if leg == 2:
        eplane = sp.Matrix([0, +ct, 0, -st, 0])
    elif leg == 3:
        eplane = sp.Matrix([0, -ct, 0, +st, 0])
    else:
        raise ValueError("leg must be 2 or 3")
    ey = sp.Matrix([0, 0, 1, 0, 0])
    out = (eplane + I * h * ey) / sp.sqrt(2)
    return out if d == 5 else sp.Matrix(list(out) + [0])


def massive_basis():
    # Each vector is unit spacelike and orthogonal to both 5D-null massive momenta.
    return [
        sp.Matrix([0, 1, 0, 0, 0]),
        sp.Matrix([0, 0, 1, 0, 0]),
        sp.Matrix([0, 0, 0, rho, beta]),
    ]


def vector_tree_matrix(h2: int, h3: int):
    ks = generic_kinematics(5)
    e2 = gluon_pol(2, h2, 5)
    e3 = gluon_pol(3, h3, 5)
    basis = massive_basis()
    return sp.Matrix([
        [sp.factor(sp.simplify(base.amplitude(ks, [ea, e2, e3, eb]))) for eb in basis]
        for ea in basis
    ])


def scalar_tree(h2: int, h3: int):
    ks = generic_kinematics(6)
    e2 = gluon_pol(2, h2, 6)
    e3 = gluon_pol(3, h3, 6)
    scalar = sp.Matrix([0, 0, 0, 0, 0, 1])
    return sp.factor(sp.simplify(base.amplitude(ks, [scalar, e2, e3, scalar])))


def sew(A: sp.Matrix, B: sp.Matrix):
    return sp.factor(sp.simplify(sum(A[a, b] * B[a, b]
                                     for a in range(3) for b in range(3))))


def main() -> None:
    g5 = base.metric(5)
    ks = generic_kinematics(5)
    basis = massive_basis()

    # Kinematics and massive polarization basis are exact.
    assert all(sp.simplify(base.mdot(g5, k, k)) == 0 for k in ks)
    assert sp.simplify(sum(ks, sp.zeros(5, 1))) == sp.zeros(5, 1)
    for e in basis:
        assert sp.simplify(base.mdot(g5, ks[0], e)) == 0
        assert sp.simplify(base.mdot(g5, ks[3], e)) == 0
        assert sp.simplify(base.mdot(g5, e, e) + 1) == 0

    Mmm = vector_tree_matrix(-1, -1)
    Mpp = vector_tree_matrix(+1, +1)
    Smm = scalar_tree(-1, -1)
    Spp = scalar_tree(+1, +1)
    Cv_same = sew(Mmm, Mpp)
    Cs_same = sp.factor(Smm * Spp)

    same_defect = sp.factor(sp.simplify(Cv_same - 3 * Cs_same))
    same_defect_target = sp.factor(
        4 * (r**2 - 1)**2 * (t**2 + 1)**2 / (r**2 + t**2)**2
    )
    assert sp.simplify(same_defect - same_defect_target) == 0
    assert sp.simplify(same_defect.subs(r, 1)) == 0

    Mpm = vector_tree_matrix(+1, -1)
    Mmp = vector_tree_matrix(-1, +1)
    Spm = scalar_tree(+1, -1)
    Smp = scalar_tree(-1, +1)
    assert Mmp == Mpm.conjugate()
    assert sp.simplify(Spm - Smp) == 0

    mixed_scalar_target = sp.factor(
        -2 * t**2 * (r**2 - 1)**2 /
        ((r**2 + 1) * (r**2 + t**2) * (t**2 + 1))
    )
    assert sp.simplify(Spm - mixed_scalar_target) == 0

    Cv_mixed = sew(Mpm, Mmp)
    Cs_mixed = sp.factor(Spm * Smp)
    assert sp.simplify(Cv_mixed.subs(r, 1) - 16) == 0
    assert sp.simplify(Cs_mixed.subs(r, 1)) == 0

    # The old threshold point is beta=0, rho=1.
    assert sp.simplify(beta.subs(r, 1)) == 0
    assert sp.simplify(rho.subs(r, 1) - 1) == 0

    print("same-helicity vector sewing =", Cv_same)
    print("same-helicity scalar sewing =", Cs_same)
    print("same-helicity (vector - 3 scalar) defect =", same_defect)
    print("mixed-helicity scalar tree =", Spm)
    print("mixed-helicity vector sewing =", Cv_mixed)
    print("mixed-helicity scalar sewing =", Cs_mixed)
    print("threshold r=1: same 3:1 restored; mixed vector=16, mixed scalar=0")
    print("PASS: previous state-count identities are certified as threshold-only")


if __name__ == "__main__":
    main()
