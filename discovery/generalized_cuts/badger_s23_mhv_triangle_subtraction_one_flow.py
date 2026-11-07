#!/usr/bin/env python3
"""Exact surviving s23 triangle subtraction for one scalar-flow orientation.

This companion uses the same S1=s23=1, Q=1 rational bispinor frame as
``badger_s23_mhv_pure_bubble_boundary.py``.  The pure mu2 bubble boundary was
shown there to vanish.  Here we compute the nontrivial triangle subtraction
associated with the uncut propagator in the right tree, K3=p1.

Key exact steps:

* D_R=(l1+p1)^2-mu2=-P/[t(1+u^2)], where P is quadratic in y.
* On P=0, the mixed four-point numerator collapses to
      N=t(1+u^2)(t+u y-u).
* The product of the two Badger three-point scalar trees (eq. 56), using
  references xi_4=p1 and xi_1=p4, obeys
      A3(q_R,4+,-k) A3(k,1-,l1) = i D_R A_R.
  Thus the genuine triple-cut product for one scalar-flow orientation is
      A_L * i D_R * A_R.
* The sum over the two y roots is performed branch-free via Vieta identities.
* Badger's T1,T2,T3 moments (eqs. 50--52) are inserted for K3=p1.

The complete moment-mapped subtraction before the Badger -1/2 prefactor reduces
exactly to

    Tmapped = i (10 mu2 u^2 + 6 mu2 - u^2) / (3(1+u^2)).

Hence its mu2 coefficient is

    [mu2] Tmapped = 2 i (5u^2+3)/(3(1+u^2)),

and after the -1/2 prefactor in Badger eq. (43), one scalar-flow orientation gives

    C_tri,one-flow^[2]
      = - i (5 u^2 + 3) / (3(1+u^2)).

The published full coefficient in this frame is twice this value.  This file
therefore isolates the remaining normalization question to the scalar-flow
multiplicity; it does NOT silently insert that factor of two.
"""

from __future__ import annotations

import sympy as sp

u, y, t, mu2 = sp.symbols("u y t mu2", nonzero=True)
I = sp.I
eps = sp.Matrix([[0, 1], [-1, 0]])

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
    return sp.factor((lama.T * eps * P * eps * tilb)[0])


def square_P_angle(tila, P, lamb):
    return sp.factor((tila.T * eps * P.T * eps * lamb)[0])


def sq(P):
    return sp.factor(P.det())


def inf_polynomial(expr, var):
    num, den = sp.fraction(sp.cancel(expr))
    q, _r = sp.div(sp.Poly(num, var), sp.Poly(den, var))
    return sp.factor(q.as_expr())


def main() -> None:
    qL = -l1
    AL = sp.factor(I * square_P_angle(til3, qL, lam2)**2 /
                   square_P_angle(til2, qL, lam2))
    assert sp.simplify(AL - I * t**2 / y) == 0

    qR = -l2
    AR = sp.factor(I * angle_P_square(lam1, qR, til4)**2 /
                   angle_P_square(lam4, qR, til4))

    DR = sp.factor(sq(l1 + p1) - mu2)
    P = sp.factor(-DR * t * (1 + u**2))
    target_P = sp.expand(
        u * y**2 + (t * (1 - u**2) - u) * y
        + u * mu2 - u * t**2 + u**2 * t
    )
    assert sp.simplify(P - target_P) == 0

    N = sp.factor(mu2*u**2 + t**2 + 2*t*u*y - t*u + u**2*y**2 - u**2*y)
    N_on_root = sp.rem(sp.Poly(N, y), sp.Poly(P, y)).as_expr()
    assert sp.simplify(N_on_root - t*(1+u**2)*(t+u*y-u)) == 0

    # Explicit factorization into Badger eq. (56) three-point trees.
    k = sp.simplify(qR + p4)
    A3_4p = sp.factor(I * square_P_angle(til4, qR, lam1) / angle(lam4, lam1))
    A3_1m = sp.factor(I * angle_P_square(lam1, k, til4) / square(til1, til4))
    A3prod = sp.factor(A3_4p * A3_1m)
    assert sp.simplify(A3prod - I * DR * AR) == 0

    triple = sp.factor(AL * A3prod)

    # On P=0, reduce the numerator and write the triple product as
    # i t^2 (t+u y-u)^2/y.  The sign follows from the explicit A3 product.
    triple_on_root = sp.factor(I * t**2 * (t + u*y - u)**2 / y)
    # Equality is modulo the root polynomial P.
    cleared = sp.factor((triple - triple_on_root) * y)
    rem = sp.rem(sp.Poly(sp.together(cleared).as_numer_denom()[0], y), sp.Poly(P, y)).as_expr()
    assert sp.simplify(rem) == 0

    # Branch-free sum over the two roots y_± using Vieta.
    B = t * (1 - u**2) - u
    C = u * mu2 - u * t**2 + u**2 * t
    sum_y = sp.factor(-B / u)
    sum_inv_y = sp.factor(-B / C)
    sum_triple = sp.factor(I * t**2 * (
        4*t*u - 4*u**2 + u**2*sum_y + (t-u)**2*sum_inv_y
    ))

    inf_t = sp.expand(inf_polynomial(sum_triple, t))
    target_inf_t = sp.expand(
        -I*mu2*t*u + I*mu2*t/u - I*mu2
        + I*t**3*u**3 + 2*I*t**3*u + I*t**3/u
        - 2*I*t**2*u**2 - 2*I*t**2 + I*t*u
    )
    assert sp.simplify(inf_t - target_inf_t) == 0

    # K3=p1 data: S1=1, S3=0, gamma_bar=1,
    # K1.K3=-1/2, Delta=1/4, <chi|K3|K1flat]=-u/(1+u^2).
    T1 = sp.factor(2*u/(1+u**2))
    T2 = sp.factor(3*u**2/(1+u**2)**2)
    T3 = sp.factor(u**3*(11 + 16*mu2)/(3*(1+u**2)**3))

    poly_t = sp.Poly(inf_t, t)
    mapped = sp.factor(
        poly_t.coeff_monomial(t) * T1
        + poly_t.coeff_monomial(t**2) * T2
        + poly_t.coeff_monomial(t**3) * T3
    )
    mapped_target = sp.factor(
        I * (10*mu2*u**2 + 6*mu2 - u**2) / (3*(1+u**2))
    )
    assert sp.simplify(mapped - mapped_target) == 0

    mu2_coeff_before_prefactor = sp.factor(sp.diff(mapped, mu2))
    mu2_coeff_target = sp.factor(2*I*(5*u**2 + 3)/(3*(1+u**2)))
    assert sp.simplify(mu2_coeff_before_prefactor - mu2_coeff_target) == 0

    one_flow = sp.factor(-sp.Rational(1, 2) * mu2_coeff_before_prefactor)
    expected_one_flow = sp.factor(-I * (5*u**2 + 3)/(3*(1+u**2)))
    assert sp.simplify(one_flow - expected_one_flow) == 0

    s12 = sp.factor(sq(p1+p2))
    published = sp.factor(sp.Rational(2, 3)*I*(2*s12 - 3))
    assert sp.simplify(2*one_flow - published) == 0

    print("D_R =", DR)
    print("P(y) =", P)
    print("A3*A3 = i D_R A_R: PASS")
    print("sum-root Inf_t =", sp.factor(inf_t))
    print("mapped T moments =", mapped)
    print("mu2 coefficient before -1/2 =", mu2_coeff_before_prefactor)
    print("one-flow C_tri^[2] =", one_flow)
    print("published C2^[2] in frame =", published)
    print("ratio published / one-flow =", sp.simplify(published/one_flow))
    print("PASS: one scalar-flow orientation gives exactly one half of published bubble coefficient")


if __name__ == "__main__":
    main()
