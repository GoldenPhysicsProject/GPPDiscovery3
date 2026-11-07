#!/usr/bin/env python3
"""Exact symbolic audit of the cyclic three-box symmetric-invariant hierarchy."""

import sympy as sp

s, t, u = sp.symbols("s t u")


def S(m, x, y):
    return sp.expand(
        sum(sp.factorial(j) * sp.factorial(m - j) * x**j * y ** (m - j)
            for j in range(m + 1))
    )


def invariant_expansion(m):
    e2 = -(s**2 + s*t + t**2)       # st+tu+us after u=-s-t
    e3 = -s*t*(s+t)                  # stu after u=-s-t
    target = sp.expand((S(m, s, t) + S(m, t, u) + S(m, u, s)).subs(u, -s-t))
    basis = [(a, b) for a in range(m // 2 + 1) for b in range(m // 3 + 1)
             if 2*a + 3*b == m]
    cs = sp.symbols(f"c0:{len(basis)}")
    ansatz = sp.expand(sum(c * e2**a * e3**b for c, (a, b) in zip(cs, basis)))
    equations = sp.Poly(target - ansatz, s, t).coeffs()
    sols = sp.solve(equations, cs, dict=True)
    assert len(sols) == 1
    reconstructed = sp.expand(ansatz.subs(sols[0]))
    assert sp.expand(target - reconstructed) == 0
    return basis, sols[0]


def main():
    expected = {
        0: {(0, 0): 3},
        1: {},
        2: {(1, 0): -7},
        3: {(0, 1): 30},
        4: {(2, 0): 88},
        5: {(1, 1): -1092},
        6: {(3, 0): -2700, (0, 2): 3924},
        7: {(2, 1): 66096},
        8: {(4, 0): 153216, (1, 2): -604656},
        9: {(3, 1): -6209280, (0, 3): 2043360},
    }

    # m=1 has no degree-one symmetric invariant after e1=0, and the sum vanishes.
    assert sp.expand((S(1, s, t) + S(1, t, u) + S(1, u, s)).subs(u, -s-t)) == 0

    for m in range(10):
        if m == 1:
            print("m=1: 0")
            continue
        basis, solution = invariant_expansion(m)
        got = {ab: sp.Integer(solution[c]) for c, ab in zip(sp.symbols(f"c0:{len(basis)}"), basis)}
        assert got == expected[m], (m, got, expected[m])
        print(f"m={m}: {got}")

    print("All cyclic three-box invariant identities passed exactly.")


if __name__ == "__main__":
    main()
