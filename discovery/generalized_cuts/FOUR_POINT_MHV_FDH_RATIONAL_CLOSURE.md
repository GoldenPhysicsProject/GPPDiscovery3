# Four-point adjacent-MHV pure-YM rational remainder in FDH

Codex/GPT discovery track, 2026-08-25.

This note closes the coefficient algebra left open after the supersymmetric reduction and the topology collapse.

## Inputs

For the color ordering

\[
1^-\,2^-\,3^+\,4^+,
\]

Badger's four-point D-dimensional coefficient list gives

\[
C_4^{[4]}
=2i\,\frac{\langle12\rangle[43]}
          {\langle34\rangle[21]}
=2i\,Q,
\]

where

\[
Q:=\frac{\langle12\rangle[34]}
        {\langle34\rangle[12]}
=-\Xi,
\qquad
\Xi=-\frac{\langle12\rangle[34]}
          {[12]\langle34\rangle}.
\]

The only nonzero independent `mu^2` bubble coefficient is

\[
\boxed{
C_{2;23}^{[2]}
=
\frac{2i(2s_{12}-3s_{23})\langle12\rangle^2[41]}
{3\langle14\rangle\langle23\rangle^3[21]^2[32]}.
}
\]

All `mu^2` triangle coefficients vanish for this adjacent-MHV helicity sector.

For massless internal scalar lines the universal rational integral limits give

\[
I_4[\mu^4]\to-\frac16,
\qquad
I_2[\mu^2]\to-\frac{s_{23}}6.
\]

Equivalently, specializing the cyclic rational-coefficient formula to four points, the complementary presentations `23` and `41` describe the same bubble channel, so the two `-1/12` cyclic contributions combine to `-1/6`.

Therefore

\[
\boxed{
R_4^{\rm FDH}(--++)
=-\frac16 C_4^{[4]}
 -\frac{s_{23}}6 C_{2;23}^{[2]}.
}
\]

## Four-point spinor reduction

Let

\[
s=s_{12},\qquad t=s_{23}.
\]

Momentum conservation gives the two useful identities

\[
\langle34\rangle[41]=\langle23\rangle[21],
\]

and

\[
\langle23\rangle[34]=\langle12\rangle[14].
\]

Together with

\[
t=\langle23\rangle[32]
=\langle14\rangle[41],
\]

these imply

\[
\boxed{
 t\,C_{2;23}^{[2]}
 = iQ\,\frac{2(2s-3t)}{3t}.
}
\]

Substitution into the rational remainder yields a complete cancellation of the `3t` pieces:

\[
\begin{aligned}
R_4^{\rm FDH}(--++)
&=-\frac{iQ}{3}
  -\frac{iQ}{9t}(2s-3t)\\
&=\boxed{-\frac{2i}{9}\frac{s}{t}\,Q}.
\end{aligned}
\]

Using `Xi=-Q`, the convention-fixed celestial-helicity form is

\[
\boxed{
R_4^{\rm FDH}(1^-,2^-,3^+,4^+)
=
\frac{2i}{9}\frac{s_{12}}{s_{23}}\,\Xi
=
-\frac{2i}{9}\frac{s_{12}}{s_{23}}
\frac{\langle12\rangle[34]}
     {[12]\langle34\rangle}.
}
\]

This is the rational contribution carried by the complex-adjoint-scalar sector and, by the supersymmetric decomposition, the pure-gluon rational remainder in the FDH organization.

## Celestial radial decomposition

The integrated rational number is assembled from two distinct hyperbolic shells before dimension shifting:

- box: `mu^4`, radial shape `tanh r sech^4 r`, peak `mu=M/sqrt(5)`;
- bubble: `mu^2`, radial shape `tanh r sech^2 r`, peak `mu=M/sqrt(6)`.

Thus the compact rational answer hides two geometrically different evanescent shell contributions whose kinematic coefficients cancel nontrivially in the final four-point expression.

## Scheme boundary

This is an FDH result. Badger states that conversion of a gluon amplitude to the 't Hooft--Veltman scheme requires subtracting the standard finite `c_Gamma A_tree/3` scheme shift. That normalization is not folded into the boxed formula above.

Primary source: S. D. Badger, *Direct Extraction Of One Loop Rational Terms*, JHEP 01 (2009) 049, arXiv:0806.4600, especially eqs. (3.14)--(3.15), (5.1)--(5.3), (5.10), and (5.33).
