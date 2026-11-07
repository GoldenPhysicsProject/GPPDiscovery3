# Exact dominant radial shell of the D-dimensional mu^4 box sector

Codex/GPT discovery track, 2026-08-25.

## Setup

For the D-dimensional two-particle cut with channel mass `M` and transverse mass `mu`,

\[
\tanh r=\sqrt{1-\frac{4\mu^2}{M^2}},\qquad
\mu=\frac{M}{2\cosh r}.
\]

The massive-scalar rational-box sector carries a numerator `mu^4`, while the integrated two-body phase space contributes a factor proportional to `tanh r`. Therefore the universal radial weight multiplying the angular two-propagator kernel is

\[
W(r)=\tanh r\,\operatorname{sech}^4 r,
\]

up to the overall factor `M^4/16` and convention-dependent phase-space normalization.

## Exact maximum

Set

\[
x=\tanh r,\qquad 0\le x<1.
\]

Since `sech^2 r = 1-x^2`,

\[
W(x)=x(1-x^2)^2.
\]

Differentiating,

\[
W'(x)=(1-x^2)(1-5x^2).
\]

Hence the unique interior critical point is

\[
\boxed{x=\tanh r_*=\frac1{\sqrt5}}.
\]

Because `W(0)=0`, `W(x)>0` for `0<x<1`, and `W(x)\to0` as `x\to1^-`, this critical point is the unique global maximum.

At the maximum,

\[
\operatorname{sech}^2 r_*=1-\frac15=\frac45,
\qquad
\cosh r_*=\frac{\sqrt5}{2},
\]

so

\[
\boxed{\mu_*=\frac{M}{\sqrt5}}.
\]

Equivalently,

\[
\boxed{\frac{4\mu_*^2}{M^2}=\frac45}.
\]

The maximal dimensionless radial weight is

\[
\boxed{W(r_*)=\frac{16}{25\sqrt5}}.
\]

Thus the phase-space-weighted numerator itself peaks at

\[
\boxed{
\mu^4\tanh r
\;\le\;
\frac{M^4}{25\sqrt5},
}
\]

with equality at `mu=M/sqrt(5)`.

## Interpretation

This is a universal statement about the isolated `mu^4` box sector before the angular propagator kernel is included. The rational numerator is suppressed both at threshold (`r=0`, `mu=M/2`, where phase space closes) and toward the four-dimensional boundary (`r->infinity`, `mu->0`, where `mu^4` vanishes). The numerator-times-phase-space measure is therefore concentrated on a finite hyperbolic shell, with its exact maximum at `tanh r=1/sqrt(5)`.

This does **not** imply that the complete cut integral is maximized at the same shell: the remaining angular kernel `J(r)` can shift the dominant support. It does provide an exact benchmark for any subsequent radial integration or saddle analysis of the D-dimensional rational sector.

## Next boundary

Multiply this universal weight by the exact massive `S^2` two-propagator master kernel and study the resulting radial integrand, including its singular/threshold structure. Only after the complete D-dimensional state sum and subtraction sectors are assembled can this be promoted to a pure-Yang--Mills rational-amplitude statement.
