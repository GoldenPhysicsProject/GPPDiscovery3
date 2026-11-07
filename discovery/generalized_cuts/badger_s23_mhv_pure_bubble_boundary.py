#!/usr/bin/env python3
"""Exact pure-bubble boundary audit for C_{2;23}^{[2]}(--++).

This implements Badger's four-point s23 double-cut parametrization (0806.4600,
eqs. 39--45) in an explicit rational bispinor frame with S1=s23=1 and chi=p2.
The external momenta are chosen as

    p2 = |2><2| = diag(1,0),
    p3 = |3><3| = diag(0,1),
    p1 + p4 = -(p2+p3),

with a rational external parameter u.  This frame has

    s23 = 1,
    s12 = -u^2/(1+u^2),
    Q = <12>[34]/(<34>[12]) = 1.

The cut loop bispinor is Badger eq. (39):

    l1 = [[1-y, (y(1-y)-mu2)/t],
          [t,   y                    ]].

Using Badger's mixed scalar tree A4(1s,2+,3-,4s) from eq. (57) and its parity
conjugate for (2-,3+), with the cut ordering

    A_L(-l1,2-,3+,l2) A_R(-l2,4+,1-,l1),

we form the exact double-cut product.  Polynomial division in y computes Inf_y.
The result is

    Inf_y[A_L A_R]
      = -t*u^2/(1+u^2) * (t*u^2 + 3*t + u*y - u),

which contains no mu2.  Consequently the coefficient selected by
Inf_{mu^2}|_{mu^2} in Badger eq. (42) vanishes:

    C_{2;23}^{bub,[2]}(--++) = 0

in this parametrization.  The known nonzero full bubble coefficient must therefore
be supplied entirely by the triangle-subtraction terms of eqs. (43)--(44).

Overall tree/cut normalization is kept in Badger's i-convention inside this audit.
"""

from __future__ import annotations

import sympy as sp

u, y, t, mu2 = sp.symbols("u y t mu2", nonzero=True)
I = sp.I
eps = sp.Matrix([[0, 1], [-1, 0]])

# Massless spinors and momenta.  A momentum matrix is lambda * tilde_lambda^T.
lam2 = sp.Matrix([1, 0]);  til2 = sp.Matrix([1, 0])
lam3 = sp.Matrix([0, 1]);  til3 = sp.Matrix([0, 1])
lam1 = sp.Matrix([1, u]);  til1 = -sp.Matrix([1, u]) / (1 + u**2)
lam4 = sp.Matrix([u, -1]); til4 = -sp.Matrix([u, -1]) / (1 + u**2)

p2 = lam2 * til2.T
p3 = lam3 * til3.T
p1 = lam1 * til1.T
p4 = lam4 * til4.T
K1 = p2 + p3

l1 = sp.Matrix([
    [1 - y, (y * (1 - y) - mu2) / t],
    [t, y],
])
l2 = l1 - K1


def angle(a, b):
    return sp.factor((a.T * eps * b)[0])


def square(a, b):
    return sp.factor((a.T * eps * b)[0])


def angle_P_square(lama, P, tilb):
    # <a|P|b]
    return sp.factor((lama.T * eps * P * eps * tilb)[0])


def square_P_angle(tila, P, lamb):
    # [a|P|b>
    return sp.factor((tila.T * eps * P.T * eps * lamb)[0])


def minkowski_sq(P):
    return sp.factor(P.det())


def inf_polynomial(expr, var):
    """Polynomial part at var=infinity of an exact rational function."""
    num, den = sp.fraction(sp.cancel(expr))
    q, _r = sp.div(sp.Poly(num, var), sp.Poly(den, var))
    return sp.factor(q.as_expr())


def main() -> None:
    assert sp.simplify(p1 + p2 + p3 + p4) == sp.zeros(2)
    assert all(sp.simplify(minkowski_sq(p)) == 0 for p in (p1, p2, p3, p4))
    assert sp.simplify(minkowski_sq(l1) - mu2) == 0
    assert sp.simplify(minkowski_sq(l2) - mu2) == 0

    s23 = sp.factor(minkowski_sq(p2 + p3))
    s12 = sp.factor(minkowski_sq(p1 + p2))
    Q = sp.factor(angle(lam1, lam2) * square(til3, til4) /
                  (angle(lam3, lam4) * square(til1, til2)))
    assert s23 == 1
    assert sp.simplify(s12 + u**2 / (1 + u**2)) == 0
    assert sp.simplify(Q - 1) == 0

    # Left tree: parity conjugate of Badger eq. (57), A4(1s,2-,3+,4s).
    qL = -l1
    AL = sp.factor(I * square_P_angle(til3, qL, lam2)**2 /
                   square_P_angle(til2, qL, lam2))

    # Right tree: Badger eq. (57), A4(1s,2+,3-,4s).
    qR = -l2
    AR = sp.factor(I * angle_P_square(lam1, qR, til4)**2 /
                   angle_P_square(lam4, qR, til4))

    assert sp.simplify(AL - I * t**2 / y) == 0

    product = sp.factor(AL * AR)
    inf_y = inf_polynomial(product, y)
    target_inf_y = sp.factor(
        -t * u**2 * (t * u**2 + 3 * t + u * y - u) / (1 + u**2)
    )
    assert sp.simplify(inf_y - target_inf_y) == 0

    # Badger Y0=1, Y1=1/2.  There is no y^2 term here.
    after_Y = sp.factor(inf_y.coeff(y, 0) + sp.Rational(1, 2) * inf_y.coeff(y, 1))
    # The entire Inf_y polynomial boundary is mu2-independent, before or after Y moments.
    assert sp.simplify(sp.diff(inf_y, mu2)) == 0
    assert sp.simplify(sp.diff(after_Y, mu2)) == 0

    inf_t = inf_polynomial(after_Y, t)
    assert sp.simplify(sp.diff(inf_t, mu2)) == 0

    full_target = sp.factor(sp.Rational(2, 3) * I * (2 * s12 - 3 * s23))
    assert sp.simplify(full_target) != 0

    print("s12 =", s12, "s23 =", s23, "Q =", Q)
    print("A_L =", AL)
    print("A_R =", AR)
    print("Inf_y[A_L A_R] =", inf_y)
    print("after Y moments =", after_Y)
    print("Inf_t boundary =", inf_t)
    print("mu2 coefficient of pure-bubble boundary = 0")
    print("known full target in Q=1 frame =", full_target)
    print("PASS: nonzero C2[2] must arise entirely from triangle-subtraction terms")


if __name__ == "__main__":
    main()
