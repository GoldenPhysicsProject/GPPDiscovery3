# Four-point adjacent-MHV bubble coefficient: phase normal form and subtraction boundary

Codex/GPT discovery track, 2026-08-26.

## 1. Literature target, rewritten without importing the final rational amplitude

Badger's explicit adjacent-MHV bubble coefficient is

\[
C_{2;23}^{[2]}(1^-,2^-,3^+,4^+)
=
\frac{2i(2s_{12}-3s_{23})\langle12\rangle^2[41]}
{3\langle14\rangle\langle23\rangle^3[21]^2[32]}.
\]

At four massless points, momentum conservation gives

\[
\langle12\rangle[32]=-\langle14\rangle[34],
\qquad
\langle34\rangle[41]=-\langle23\rangle[12].
\]

Using also

\[
s_{23}=\langle23\rangle[32],
\qquad
Q:=\frac{\langle12\rangle[34]}{\langle34\rangle[12]},
\]

the spinor prefactor collapses exactly:

\[
\frac{\langle12\rangle^2[41]}
{\langle14\rangle\langle23\rangle^3[21]^2[32]}
=
\frac{Q}{s_{23}^2}.
\]

Therefore the bubble target has the phase-normal form

\[
\boxed{
C_{2;23}^{[2]}(--++)
=
\frac{2i}{3}\,
\frac{2s_{12}-3s_{23}}{s_{23}^2}\,Q.
}
\]

This is a useful target for a direct cut derivation because all helicity dependence is
isolated in the same dimensionless phase `Q` that appears in the `mu^4` box coefficient.

## 2. Important subtraction correction

Badger's double-cut extraction is not just the polynomial boundary of the double cut. In
his notation,

\[
C_2^{[2]}
=
C_2^{{\rm bub},[2]}
+
\sum_{\{K_3\}} C_2^{{\rm tri}(K_3),[2]}.
\]

The pure-bubble term is obtained from the large-`y`, large-`t`, large-`mu^2` boundary,
with the `y` moments

\[
Y_0=1,\qquad
Y_1=\frac12,\qquad
Y_2=\frac13\left(1-\frac{\mu^2}{S_1}\right).
\]

The triangle-subtraction pieces are evaluated from triple-cut data using nonzero powers
of the free `t` parameter and the corresponding `T_i` moments.

Consequently, the fact that the *final triangle integral coefficients*
`C_3^[2](--++)` vanish does **not** imply that the triangle-subtraction contribution to
the bubble extraction vanishes.  The former probes the `t^0` boundary of the triple cut;
the latter can depend on `t^1`, `t^2`, and higher moments.  Any direct derivation of the
`23` bubble must therefore retain and evaluate the higher-topology residues before
claiming the coefficient.

## 3. State-sum audit correction

The earlier exact `3:1` massive-vector/scalar state sum and the mixed-helicity `16` state
sum were evaluated on a threshold slice in which the four-dimensional projections of the
massive legs are at rest.  The generic rational two-parameter audit now shows that these
relations are not generic cut identities.  They remain exact threshold regression tests
only and cannot be used to normalize the bubble coefficient.

## Boundary

The boxed phase-normal form above is an algebraic rewriting of the known literature
coefficient, not yet an independent cut derivation.  The next independent step is to
construct the generic `23` double-cut product of massive trees, take its prescribed
`Inf_y`, `Inf_t`, and `Inf_mu^2` boundary, and add the required triple-cut subtraction
moments.  Only after that result matches the boxed expression is the bubble channel
independently closed.

External source checked: S. D. Badger, JHEP 01 (2009) 049, arXiv:0806.4600, especially
eqs. (4.3)--(4.10), (5.1)--(5.3), and the four-point bubble list containing the
adjacent-MHV `23` coefficient.
