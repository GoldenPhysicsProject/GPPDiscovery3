#!/usr/bin/env python3
"""Exact symbolic adjacent-MHV massive-vector state sum at threshold.

Codex/GPT Golden Physics discovery track, 2026-08-26.

This script uses the rational angular parametrization

    cos(theta) = (1-t^2)/(1+t^2),
    sin(theta) = 2t/(1+t^2),

and keeps the centre-of-mass energy E symbolic.  It exactly evaluates the complete
color-ordered four-gluon tree (two cubic channels plus quartic contact term) when the
5D-null internal legs have four-dimensional projections at rest.  On this threshold
slice the same-helicity tree matrices reduce to

    A_h^{ab} = -delta^{ab},   h=+/- , a,b=1,2,3,

while the extra-dimensional scalar tree is

    A_h^(S) = +1.

Consequently, at threshold only,

    C^(V_m) = 3 C^(S),
    C^(V_m) - C^(S) = 2 C^(S).

Important correction: ``massive_vector_generic_state_sum_symbolic.py`` keeps an
independent rational mass/velocity parameter and proves that the 3:1 relation is not a
generic cut identity.  This file is retained as an exact threshold regression test, not
as a derivation of a generic Yang--Mills cut coefficient.
"""

from __future__ import annotations

import sympy as sp

E, t = sp.symbols("E t", real=True, nonzero=True)
I = sp.I
c = (1 - t**2) / (1 + t**2)
s = 2 * t / (1 + t**2)


def metric(d: int) -> sp.Matrix:
    return sp.diag(1, *([-1] * (d - 1)))


def mdot(g: sp.Matrix, a: sp.Matrix, b: sp.Matrix):
    return (a.T * g * b)[0]


def v3(g: sp.Matrix, p: sp.Matrix, q: sp.Matrix, r: sp.Matrix):
    d = len(p)
    out = sp.MutableDenseNDimArray.zeros(d, d, d)
    for mu in range(d):
        for nu in range(d):
            for rho in range(d):
                out[mu, nu, rho] = (
                    g[mu, nu] * (p - q)[rho]
                    + g[nu, rho] * (q - r)[mu]
                    + g[rho, mu] * (r - p)[nu]
                )
    return out


def contract_12_34(e1, e2, V1, g, V2, e3, e4):
    d = len(e1)
    total = 0
    for m in range(d):
        for n in range(d):
            if e1[m] == 0 or e2[n] == 0:
                continue
            for a in range(d):
                for b in range(d):
                    if g[a, b] == 0:
                        continue
                    for r in range(d):
                        if e3[r] == 0:
                            continue
                        for s_ in range(d):
                            if e4[s_] == 0:
                                continue
                            total += (
                                e1[m] * e2[n] * V1[m, n, a] * g[a, b]
                                * V2[b, r, s_] * e3[r] * e4[s_]
                            )
    return sp.factor(total)


def contract_23_41(e2, e3, V1, g, V2, e4, e1):
    d = len(e1)
    total = 0
    for n in range(d):
        for r in range(d):
            if e2[n] == 0 or e3[r] == 0:
                continue
            for a in range(d):
                for b in range(d):
                    if g[a, b] == 0:
                        continue
                    for s_ in range(d):
                        if e4[s_] == 0:
                            continue
                        for m in range(d):
                            if e1[m] == 0:
                                continue
                            total += (
                                e2[n] * e3[r] * V1[n, r, a] * g[a, b]
                                * V2[b, s_, m] * e4[s_] * e1[m]
                            )
    return sp.factor(total)


def amplitude(ks: list[sp.Matrix], eps: list[sp.Matrix]):
    d = len(ks[0])
    g = metric(d)
    e1, e2, e3, e4 = [g * e for e in eps]
    k1, k2, k3, k4 = ks

    p12 = k1 + k2
    a12 = contract_12_34(
        e1, e2, v3(g, k1, k2, -p12), g,
        v3(g, p12, k3, k4), e3, e4,
    ) / mdot(g, p12, p12)

    p23 = k2 + k3
    a23 = contract_23_41(
        e2, e3, v3(g, k2, k3, -p23), g,
        v3(g, p23, k4, k1), e4, e1,
    ) / mdot(g, p23, p23)

    c12_34 = mdot(g, eps[0], eps[1]) * mdot(g, eps[2], eps[3])
    c13_24 = mdot(g, eps[0], eps[2]) * mdot(g, eps[1], eps[3])
    c14_23 = mdot(g, eps[0], eps[3]) * mdot(g, eps[1], eps[2])
    return sp.factor(a12 + a23 - c12_34 + 2 * c13_24 - c14_23)


def kinematics(d: int = 5) -> list[sp.Matrix]:
    # The four-dimensional projections of legs 1 and 2 are both (-E,0,0,0):
    # this is the threshold/rest slice.  Their opposite fifth components supply mass E.
    base = [
        sp.Matrix([-E, 0, 0, 0, -E]),
        sp.Matrix([-E, 0, 0, 0, +E]),
        sp.Matrix([+E, E * c, E * s, 0, 0]),
        sp.Matrix([+E, -E * c, -E * s, 0, 0]),
    ]
    if d == 5:
        return base
    if d == 6:
        return [sp.Matrix(list(k) + [0]) for k in base]
    raise ValueError("only d=5 or d=6 is used")


def gluon_helicity(leg: int, h: int, d: int = 5) -> sp.Matrix:
    if leg == 3:
        eplane = sp.Matrix([0, -s, c, 0, 0])
    elif leg == 4:
        eplane = sp.Matrix([0, +s, -c, 0, 0])
    else:
        raise ValueError("leg must be 3 or 4")
    ez = sp.Matrix([0, 0, 0, 1, 0])
    out = (eplane + I * h * ez) / sp.sqrt(2)
    return out if d == 5 else sp.Matrix(list(out) + [0])


def vector_tree_matrix(h: int) -> sp.Matrix:
    basis = [
        sp.Matrix([0, 1, 0, 0, 0]),
        sp.Matrix([0, 0, 1, 0, 0]),
        sp.Matrix([0, 0, 0, 1, 0]),
    ]
    ks = kinematics(5)
    e3 = gluon_helicity(3, h, 5)
    e4 = gluon_helicity(4, h, 5)
    return sp.Matrix([
        [sp.factor(amplitude(ks, [ea, eb, e3, e4])) for eb in basis]
        for ea in basis
    ])


def scalar_tree(h: int):
    ks = kinematics(6)
    e3 = gluon_helicity(3, h, 6)
    e4 = gluon_helicity(4, h, 6)
    scalar = sp.Matrix([0, 0, 0, 0, 0, 1])
    return sp.factor(amplitude(ks, [scalar, scalar, e3, e4]))


def main() -> None:
    target = -sp.eye(3)
    matrices = {h: vector_tree_matrix(h) for h in (-1, +1)}
    scalars = {h: scalar_tree(h) for h in (-1, +1)}

    for h in (-1, +1):
        assert matrices[h] == target, (h, matrices[h])
        assert scalars[h] == 1, (h, scalars[h])

    cv = sp.expand(sum(matrices[-1][a, b] * matrices[+1][a, b]
                       for a in range(3) for b in range(3)))
    cs = sp.expand(scalars[-1] * scalars[+1])
    assert sp.simplify(cv - 3 * cs) == 0

    print("THRESHOLD ONLY: four-dimensional massive legs are at rest")
    print("A_- =")
    sp.pprint(matrices[-1])
    print("A_+ =")
    sp.pprint(matrices[+1])
    print(f"A_-^(S) = {scalars[-1]}")
    print(f"A_+^(S) = {scalars[+1]}")
    print(f"C^(V_m) = {cv}, C^(S) = {cs}")
    print("PASS: exact threshold identity C^(V_m) = 3 C^(S).")
    print("Generic cut identity is tested separately and is not 3:1.")


if __name__ == "__main__":
    main()
