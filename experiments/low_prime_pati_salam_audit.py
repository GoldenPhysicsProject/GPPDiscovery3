#!/usr/bin/env python3
"""Exact audit of the low-prime 2-3-5-7 Pati-Salam candidate.

The calculation uses exact SymPy arithmetic.  It verifies:
  * the minimal-gap path 2--3--5--7 has A3 Cartan incidence Gram;
  * its gap-2 subpath 3--5--7 has A2 Cartan incidence Gram;
  * nearest-neighbour matrix units and adjoints Lie-generate sl(4) (dim 15),
    while the gap-2 color subset Lie-generates sl(3) (dim 8);
  * the primitive trace-zero diagonal constant on 3,5,7 is X=(-3,1,1,1);
  * X/3 gives B-L=(-1,1/3,1/3,1/3);
  * Pati-Salam Q=T3L+T3R+X/6 reproduces one-generation charges.

This is an arithmetic/representation-theory audit, not a derivation of the
physical identification.
"""

from __future__ import annotations

from itertools import combinations
import sympy as sp


PRIMES = (2, 3, 5, 7)


def incidence(n: int) -> sp.Matrix:
    """Oriented incidence of an n-vertex path."""
    B = sp.zeros(n, n - 1)
    for j in range(n - 1):
        B[j, j] = 1
        B[j + 1, j] = -1
    return B


def E(n: int, i: int, j: int) -> sp.Matrix:
    M = sp.zeros(n)
    M[i, j] = 1
    return M


def flatten(M: sp.Matrix) -> sp.Matrix:
    return sp.Matrix(list(M))


def span_rank(mats: list[sp.Matrix]) -> int:
    if not mats:
        return 0
    return sp.Matrix.hstack(*(flatten(M) for M in mats)).rank()


def add_if_independent(basis: list[sp.Matrix], M: sp.Matrix) -> bool:
    if M == sp.zeros(*M.shape):
        return False
    r0 = span_rank(basis)
    r1 = span_rank(basis + [M])
    if r1 > r0:
        basis.append(M)
        return True
    return False


def lie_closure(gens: list[sp.Matrix]) -> list[sp.Matrix]:
    basis: list[sp.Matrix] = []
    for G in gens:
        add_if_independent(basis, G)

    changed = True
    while changed:
        changed = False
        snapshot = list(basis)
        for A, B in combinations(snapshot, 2):
            C = A * B - B * A
            if add_if_independent(basis, C):
                changed = True
    return basis


def main() -> None:
    B4 = incidence(4)
    A3 = B4.T * B4
    expected_A3 = sp.Matrix([[2, -1, 0], [-1, 2, -1], [0, -1, 2]])

    B3 = incidence(3)
    A2 = B3.T * B3
    expected_A2 = sp.Matrix([[2, -1], [-1, 2]])

    print("low primes:", PRIMES)
    print("successive gaps:", tuple(PRIMES[i + 1] - PRIMES[i] for i in range(3)))
    print("A3 incidence Gram =")
    print(A3)
    print("A3 exact:", A3 == expected_A3, "det =", A3.det())
    print("A2 gap-2 incidence Gram =")
    print(A2)
    print("A2 exact:", A2 == expected_A2, "det =", A2.det())
    print()

    # Ordered basis: 0=lepton candidate p=2, 1..3=color candidates p=3,5,7.
    color_gens = [E(4, 1, 2), E(4, 2, 1), E(4, 2, 3), E(4, 3, 2)]
    full_gens = [E(4, 0, 1), E(4, 1, 0)] + color_gens
    color_alg = lie_closure(color_gens)
    full_alg = lie_closure(full_gens)
    print("Lie closure dimension from gap-2 color simple roots:", len(color_alg))
    print("Lie closure dimension after adding gap-1 root:", len(full_alg))
    assert len(color_alg) == 8
    assert len(full_alg) == 15

    E35, E57 = E(4, 1, 2), E(4, 2, 3)
    E37 = E35 * E57 - E57 * E35
    print("[E_35,E_57] = E_37:", E37 == E(4, 1, 3))

    E23 = E(4, 0, 1)
    E25 = E23 * E35 - E35 * E23
    E27 = E25 * E57 - E57 * E25
    print("[E_23,E_35] = E_25:", E25 == E(4, 0, 2))
    print("[E_25,E_57] = E_27:", E27 == E(4, 0, 3))
    print()

    X = sp.Matrix([-3, 1, 1, 1])
    print("X = 3(B-L):", tuple(X))
    print("trace X =", sum(X))
    print("color-root X differences:", X[1] - X[2], X[2] - X[3])
    print("lepton-color X difference:", X[0] - X[1])
    print("B-L:", tuple(sp.Rational(v, 3) for v in X))
    print("adjoint branching dimensions: 15 = 8 + 1 + 3 + 3")
    print()
    print("SU(4) root charges under X and B-L")
    for i, j, label in [
        (1, 2, "3<->5 color"),
        (2, 3, "5<->7 color"),
        (1, 3, "3<->7 composite color"),
        (0, 1, "2<->3 leptoquark"),
        (0, 2, "2<->5 leptoquark"),
        (0, 3, "2<->7 leptoquark"),
    ]:
        dx = X[i] - X[j]
        print(f"{label:24s}: Delta X={dx:2d}, Delta(B-L)={sp.Rational(dx,3)}")
    print("Thus color roots are B-L neutral; lepton-color roots carry |B-L|=4/3.")
    print()


    # One-generation Pati-Salam charge table.
    def Q(x: int, tL: sp.Rational, tR: sp.Rational) -> sp.Rational:
        return sp.simplify(tL + tR + sp.Rational(x, 6))

    rows = [
        ("nu_L", -3, sp.Rational(1, 2), 0),
        ("e_L", -3, sp.Rational(-1, 2), 0),
        ("u_L", 1, sp.Rational(1, 2), 0),
        ("d_L", 1, sp.Rational(-1, 2), 0),
        ("nu_R", -3, 0, sp.Rational(1, 2)),
        ("e_R", -3, 0, sp.Rational(-1, 2)),
        ("u_R", 1, 0, sp.Rational(1, 2)),
        ("d_R", 1, 0, sp.Rational(-1, 2)),
    ]
    print("one-generation charge table")
    for name, x, tL, tR in rows:
        y = sp.simplify(tR + sp.Rational(x, 6))
        print(f"{name:4s}: X={x:2d}, T3L={str(tL):>4s}, T3R={str(tR):>4s}, "
              f"Y={str(y):>4s}, Q={str(Q(x,tL,tR)):>4s}")

    # Minimal-gap uniqueness facts.
    # Among primes, gap 1 can occur only at (2,3); gap 2 triple can only be 3,5,7.
    print()
    print("Interpretive checkpoint:")
    print("  p=2 is singled out as the SU(3)-singlet coordinate.")
    print("  p=3,5,7 form the A2 color candidate.")
    print("  the unique 2-3 root extends A2 to A3.")
    print("  This is a candidate dictionary; the matrix identities above are exact.")


if __name__ == "__main__":
    main()
