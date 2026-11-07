#!/usr/bin/env python3
"""Exact symbolic mixed-helicity massive-vector threshold audit.

This companion imports the convention-locked tree engine from
``massive_vector_mhv_state_sum_symbolic.py`` and evaluates the two massive-vector / two
gluon tree for opposite external helicities on the same special kinematic slice used by
that original audit.

Important correction: in this slice the four-dimensional projections of the two massive
legs are at rest.  Equivalently, this is the ``r=1`` threshold point of the generic
parametrization in ``massive_vector_generic_state_sum_symbolic.py``.  Therefore the
identities proved here are exact threshold identities, not generic two-particle-cut
identities and not bubble coefficients.

At threshold:

* the extra-dimensional adjoint-scalar mixed-helicity tree vanishes;
* the massive-vector tree matrix is rank one and factorizes as ``2 w w^T``;
* ``w^T w = 0`` and ``w^† w = 2``;
* sewing the conjugate mixed-helicity matrices gives the constant 16.

The generic audit shows that away from threshold the mixed scalar tree is nonzero and the
vector sewing is kinematics dependent.  This file is retained as a useful exact
cross-check of the threshold limit only.
"""

from __future__ import annotations

import sympy as sp

import massive_vector_mhv_state_sum_symbolic as base


def vector_tree_matrix(h3: int, h4: int) -> sp.Matrix:
    basis = [
        sp.Matrix([0, 1, 0, 0, 0]),
        sp.Matrix([0, 0, 1, 0, 0]),
        sp.Matrix([0, 0, 0, 1, 0]),
    ]
    ks = base.kinematics(5)
    e3 = base.gluon_helicity(3, h3, 5)
    e4 = base.gluon_helicity(4, h4, 5)
    return sp.Matrix([
        [sp.factor(base.amplitude(ks, [ea, eb, e3, e4])) for eb in basis]
        for ea in basis
    ])


def scalar_tree(h3: int, h4: int):
    ks = base.kinematics(6)
    e3 = base.gluon_helicity(3, h3, 6)
    e4 = base.gluon_helicity(4, h4, 6)
    scalar = sp.Matrix([0, 0, 0, 0, 0, 1])
    return sp.factor(base.amplitude(ks, [scalar, scalar, e3, e4]))


def frobenius_sew(A: sp.Matrix, B: sp.Matrix):
    return sp.factor(sum(A[a, b] * B[a, b] for a in range(3) for b in range(3)))


def simplify_matrix(M: sp.Matrix) -> sp.Matrix:
    return M.applyfunc(lambda z: sp.factor(sp.simplify(z)))


def main() -> None:
    Mmp = vector_tree_matrix(-1, +1)
    Mpm = vector_tree_matrix(+1, -1)
    Smp = scalar_tree(-1, +1)
    Spm = scalar_tree(+1, -1)

    assert sp.simplify(Smp) == 0
    assert sp.simplify(Spm) == 0
    assert Mpm == Mmp.conjugate()

    den = 1 + base.t**2
    w = sp.Matrix([
        -2 * sp.I * base.t / den,
        -sp.I * (base.t**2 - 1) / den,
        1,
    ])
    assert simplify_matrix(Mmp - 2 * w * w.T) == sp.zeros(3)
    assert sp.simplify((w.T * w)[0]) == 0
    assert sp.simplify((w.conjugate().T * w)[0] - 2) == 0

    Cv = frobenius_sew(Mmp, Mpm)
    Cs = sp.factor(Smp * Spm)
    assert sp.simplify(Cv - 16) == 0
    assert sp.simplify(Cs) == 0

    CvT = sp.factor(sum(Mmp[a, b] * Mpm[b, a]
                        for a in range(3) for b in range(3)))
    assert sp.simplify(CvT - 16) == 0
    assert sp.simplify(4 * (w.conjugate().T * w)[0] ** 2 - 16) == 0

    print("THRESHOLD ONLY: four-dimensional massive legs are at rest")
    print("mixed-helicity scalar trees (-+),(+-) =", Smp, Spm)
    print("rank-one factorization M_-+ = 2 w w^T: PASS")
    print("w^T w =", sp.simplify((w.T * w)[0]))
    print("w^dagger w =", sp.simplify((w.conjugate().T * w)[0]))
    print("mixed-helicity vector sewing at threshold =", Cv)
    print("PASS: exact threshold identity; no generic bubble claim")


if __name__ == "__main__":
    main()
