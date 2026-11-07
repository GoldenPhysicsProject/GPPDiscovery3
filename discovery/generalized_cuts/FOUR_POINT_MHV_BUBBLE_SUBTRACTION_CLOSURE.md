# Four-point adjacent-MHV mu^2 bubble: triangle-subtraction closure

Codex/GPT discovery track, 2026-08-26.

## Result

For the scalar-loop rational sector of the color-ordered four-gluon amplitude

\[
1^-\,2^-\,3^+\,4^+,
\]

the `s23` bubble coefficient in Badger's D-dimensional/Forde organization is

\[
\boxed{
C_{2;23}^{[2]}(--++)
=\frac{2i}{3}\frac{2s_{12}-3s_{23}}{s_{23}^{2}}\,
Q,
\qquad
Q=\frac{\langle12\rangle[34]}{\langle34\rangle[12]}.
}
\]

The new point established by the companion exact symbolic audits is how this coefficient is
assembled from the double cut.

## 1. The pure-bubble mu^2 boundary vanishes

Use Badger's recommended four-point choice

\[
K_1=p_2+p_3,\qquad \chi=p_2.
\]

In the rational `s23=1`, `Q=1` bispinor frame of
`badger_s23_mhv_pure_bubble_boundary.py`, the exact double-cut product has

\[
\operatorname{Inf}_y[A_LA_R]
=-\frac{t u^2}{1+u^2}
\left(tu^2+3t+uy-u\right).
\]

This polynomial is independent of `mu2`.  It remains `mu2`-independent after the `Y_i`
moments and the subsequent `Inf_t` operation.  Therefore Badger eq. (42) gives

\[
\boxed{C_{2;23}^{{\rm bub},[2]}(--++)=0.}
\]

Thus the known nonzero coefficient must arise from the triangle-subtraction terms of eq. (43).

## 2. The nontrivial triangle pole is quadratic and branch-free after summation

The right-tree uncut propagator is proportional to

\[
P(y)=u y^2+\bigl[t(1-u^2)-u\bigr]y
+u\mu^2-u t^2+u^2t.
\]

Its discriminant is exactly

\[
\boxed{
\Delta_y=\bigl[t(1+u^2)-u\bigr]^2-4\mu^2u^2.
}
\]

The two Badger subtraction solutions are the roots of this one quadratic.  The exact Vieta
relations allow the sum over the two roots to be performed before any large-`t` expansion, so no
spurious square root survives.  Their leading branches are `y_+ ~ u t` and `y_- ~ -t/u`.

The companion audit `badger_s23_triangle_pole_quadratic.py` certifies these identities.

## 3. Exact three-point factorization fixes the triple-cut numerator

For the surviving channel `K3=p1`, the right mixed-helicity four-point tree factorizes at its
uncut propagator pole.  Using Badger's three-point scalar-gluon trees (eq. 56) with explicit
reference spinors, the symbolic audit proves

\[
A_3(q_R,4^+,-k)\,A_3(k,1^-,\ell_1)
=i\,D_R\,A_R,
\]

where

\[
D_R=(\ell_1+p_1)^2-\mu^2.
\]

Thus the genuine triple-cut product is obtained directly from the tree factorization, rather than
by guessing a residue normalization.

On the two roots, the numerator simplifies exactly to

\[
N=t(1+u^2)(t+uy-u).
\]

After summing the two roots branch-free, expanding at `t=infinity`, and inserting Badger's exact
`T_1,T_2,T_3` moments, a single scalar-flow orientation gives

\[
\boxed{
C_{2;23}^{\mathrm{tri},[2]}\big|_{\rm one\ flow}
=-\frac{i(5u^2+3)}{3(1+u^2)}.
}
\]

This is certified by `badger_s23_mhv_triangle_subtraction_one_flow.py`.

## 4. Complex-scalar multiplicity closes the coefficient

The `A^[s]`/`N=0` sector used in the standard supersymmetric decomposition is the contribution of
a **complex scalar**, not one real scalar.  For purely gluonic external states the two real
components, equivalently the two scalar-flow orientations, contribute equally.  Hence the
one-flow result is doubled:

\[
C_{2;23}^{[2]}
=2\,C_{2;23}^{\mathrm{tri},[2]}\big|_{\rm one\ flow}
=-\frac{2i(5u^2+3)}{3(1+u^2)}.
\]

In the chosen rational frame,

\[
s_{23}=1,
\qquad
s_{12}=-\frac{u^2}{1+u^2},
\qquad
Q=1,
\]

so

\[
\frac{2i}{3}(2s_{12}-3)
=-\frac{2i(5u^2+3)}{3(1+u^2)},
\]

exactly reproducing the invariant target.  Restoring homogeneous dimensions and the helicity phase
therefore gives

\[
\boxed{
C_{2;23}^{[2]}(--++)
=\frac{2i}{3}\frac{2s_{12}-3s_{23}}{s_{23}^{2}}Q.
}
\]

## Structural lesson

The adjacent-MHV four-point example demonstrates explicitly that

\[
C_3^{[2]}=0
\]

does **not** imply that triangle information can be discarded from the bubble extraction.  Here the
pure bubble `mu2` boundary is zero, while the full bubble coefficient is nonzero and is supplied by
the higher-`t` moments of the triple-cut subtraction.

## Scope

This closes the massive-scalar `s23` bubble coefficient in Badger's scalar-loop rational sector and
its stated spinor/integral conventions.  It is not by itself a derivation of the complete pure-Yang--Mills
state sum, FDH-to-HV conversion, color normalization, or loop-measure normalization outside those
conventions.

Primary reference: S. D. Badger, *Direct Extraction Of One Loop Rational Terms*, arXiv:0806.4600,
especially eqs. (42)--(52), (56)--(57), and (85).  The scalar sector is the standard complex-scalar
component of the supersymmetric decomposition.
