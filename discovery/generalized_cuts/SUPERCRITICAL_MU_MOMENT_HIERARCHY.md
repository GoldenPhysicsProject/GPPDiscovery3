# Supercritical mu-moment hierarchy: from simplex volume to Symanzik moments

Codex/GPT discovery track, 2026-08-25.

The critical family previously isolated,

\[
r=n-2,
\]

is only the `m=0` member of a larger dimension-shift hierarchy. Write

\[
\boxed{r=n-2+m},\qquad m\in\mathbb Z_{\ge0}.
\]

Then the insertion is `mu^(2r)` and the raised dimension is

\[
D'=4+2r-2\epsilon
=2n+2m-2\epsilon.
\]

## Dimension-shift zero

For integer `r>=1`,

\[
\frac{\Gamma(r-\epsilon)}{\Gamma(-\epsilon)}
=(-\epsilon)(1-\epsilon)\cdots(r-1-\epsilon)
=-(r-1)!\,\epsilon+O(\epsilon^2).
\]

## Raised-dimensional ultraviolet pole

The Feynman-parameter Gamma factor is

\[
\Gamma\!\left(n-\frac{D'}2\right)
=\Gamma(-m+\epsilon).
\]

At a nonnegative integer `m`,

\[
\Gamma(-m+\epsilon)
=\frac{(-1)^m}{m!\,\epsilon}+O(1).
\]

Therefore the zero and pole combine to the universal finite coefficient

\[
\boxed{
(-1)^{m+1}\frac{(r-1)!}{m!}.
}
\]

## Parameter-space content

At `epsilon=0`, the Feynman polynomial enters with exponent

\[
\frac{D'}2-n=m.
\]

Thus, in the same scalar-integral/Feynman-polynomial convention as the critical-family derivation,

\[
\boxed{
I_n^{4-2\epsilon}[\mu^{2r}]
\longrightarrow
(-1)^{m+1}\frac{(r-1)!}{m!}
\int_{\Delta_{n-1}} F(a)^m\,da,
\qquad r=n-2+m.
}
\]

The sign of odd-`m` kinematic polynomials of course follows the chosen Minkowski/Euclidean convention for `F`; the coefficient above is the Gamma-function coefficient before any redefinition `F -> -F`.

## Critical family recovered

For `m=0`, `r=n-2`, the moment is just the simplex volume,

\[
\int_{\Delta_{n-1}}1=\frac1{(n-1)!},
\]

so

\[
I_n[\mu^{2(n-2)}]
\to
-\frac{(n-3)!}{(n-1)!}
=-\frac1{(n-1)(n-2)},
\]

exactly the earlier critical n-gon law.

## Four-graviton all-plus box

For the gravity box,

\[
n=4,\qquad r=4,\qquad m=2.
\]

The universal coefficient is

\[
(-1)^3\frac{3!}{2!}=-3.
\]

The surviving simplex datum is the quadratic moment

\[
\int_{\Delta_3}F^2
=\frac{2s^2+2t^2+st}{2520},
\]

hence

\[
I_4[\mu^8]
\to
-\frac{2s^2+2t^2+st}{840},
\]

reproducing Bern--Dixon--Perelstein--Rozowsky Appendix D.

## Interpretation

The numerator power measures how deeply the D-dimensional rational sector probes Feynman-parameter geometry:

- `m=0`: simplex volume only;
- `m=1`: first Symanzik moment;
- `m=2`: quadratic Symanzik moment;
- in general `m`: the `m`th polynomial moment of the simplex.

The four-dimensional rational residue is therefore not merely an evanescent remnant: its kinematic complexity is graded by the excess transverse numerator power above the critical value `mu^(2(n-2))`.
