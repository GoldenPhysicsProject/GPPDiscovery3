# Celestial radial shell of the mu^2 bubble rational sector

Codex/GPT discovery track, 2026-08-25.

For the four-point adjacent-MHV rational sector, the explicit D-dimensional coefficient analysis leaves only

- the already-controlled `mu^4` box sector, and
- one `mu^2` bubble sector,

with all `mu^2` triangle coefficients vanishing.

A genuine scalar bubble basis term has only the two propagators placed on the two-particle cut. After the higher-topology subtraction that defines the bubble coefficient, there are no remaining uncut propagators carrying angular dependence. Thus at fixed transverse mass `mu` the angular integration is simply

\[
\int_{S^2} d\Omega = 4\pi.
\]

With channel mass `M`, write

\[
\mu=\frac{M}{2\cosh r},
\qquad
\beta=\sqrt{1-\frac{4\mu^2}{M^2}}=\tanh r.
\]

The two-body phase-space factor contributes one power of `beta`, while the rational bubble numerator contributes `mu^2`. Therefore the universal radial shape is

\[
\boxed{
W_{\rm bub}(r)
=\tanh r\,\operatorname{sech}^2 r,
}
\]

up to the topology/helicity coefficient and the overall factor `M^2/4`.

This is exactly the `k=1` member of the universal `mu^(2k)` shell family.

## Exact normalization

Since

\[
\frac{d}{dr}\operatorname{sech}^2 r
=-2\tanh r\,\operatorname{sech}^2 r,
\]

\[
\boxed{
\int_0^\infty W_{\rm bub}(r)\,dr=\frac12.
}
\]

Hence the normalized bubble shell is

\[
\boxed{
\rho_{\rm bub}(r)
=2\tanh r\,\operatorname{sech}^2 r.
}
\]

Its CDF is

\[
F(R)=1-\operatorname{sech}^2 R=\tanh^2 R.
\]

Equivalently,

\[
\boxed{
U=\operatorname{sech}^2 r
=\left(\frac{2\mu}{M}\right)^2
\sim {\rm Uniform}(0,1).
}
\]

Thus the induced transverse-mass density is linear:

\[
\boxed{
\rho_\mu(\mu)=\frac{8\mu}{M^2},
\qquad 0\le\mu\le\frac M2.
}
\]

and

\[
\boxed{
\mathbb E[\mu^q]
=\frac{2}{q+2}\left(\frac M2\right)^q,
\qquad q>-2.
}
\]

In particular,

\[
\mathbb E\!\left[\left(\frac{2\mu}{M}\right)^2\right]=\frac12.
\]

## Unique shell maximum

For `x=tanh r`,

\[
W_{\rm bub}=x(1-x^2).
\]

Hence

\[
\frac{dW_{\rm bub}}{dx}=1-3x^2,
\]

so the unique global maximum occurs at

\[
\boxed{
\tanh r_* = \frac1{\sqrt3},
\qquad
\operatorname{sech}^2 r_* = \frac23,
\qquad
\mu_* = \frac{M}{\sqrt6}.
}
\]

The maximum value is

\[
\boxed{W_{\rm bub}^{\max}=\frac{2}{3\sqrt3}.}
\]

## Threshold and boundary behavior

At threshold `r -> 0`,

\[
W_{\rm bub}(r)=r+O(r^3),
\]

while toward the massless boundary `r -> infinity`,

\[
W_{\rm bub}(r)=4e^{-2r}+O(e^{-4r}).
\]

This is parametrically broader than the `mu^4` box shell, whose tail is `16e^{-4r}`.

## Boundary

This note closes the universal celestial geometry of the rational bubble basis insertion. It does not yet compute the adjacent-MHV coefficient `C_2^[2]` itself from the massive scalar double cut; that coefficient is the remaining helicity/kinematic datum needed before assembling the four-point rational remainder.
