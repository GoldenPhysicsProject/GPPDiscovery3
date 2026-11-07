# Four-graviton all-plus one-loop amplitude: exact mu^8 box closure

Codex/GPT discovery track, 2026-08-25.

Primary source: Bern, Dixon, Perelstein and Rozowsky, *Multi-Leg One-Loop Gravity Amplitudes from Gauge Theory*, Nucl. Phys. B546 (1999) 423, arXiv:hep-th/9811140, sec. 4.3 and appendix D.

## Scalar-loop reduction is exact in this helicity sector

For the four-graviton all-plus amplitude, the supersymmetry Ward identities permit the graviton loop to be replaced by two real minimally coupled scalars. In the two-particle `s12` cut the factor of two from the two real scalars is cancelled by the identical-particle phase-space factor `1/2`.

The required gauge-theory scalar trees are proportional to `mu^2`. Applying the four-point KLT relation squares this factor, so the corresponding gravity scalar tree is

\[
M_4^{\rm tree}(-L_1^s,1^+,2^+,L_3^s)
=-i\left(\mu^2\frac{[12]}{\langle12\rangle}\right)^2
\left[
\frac1{(\ell_1-k_1)^2-\mu^2}
+
\frac1{(\ell_1-k_2)^2-\mu^2}
\right].
\]

Thus each side of the cut carries `mu^4`, and their product carries

\[
\boxed{\mu^8}.
\]

Bern et al. write the complete `s12` cut explicitly with this `mu^8` numerator (their eq. 4.11).

## Complete four-point integral representation

Combining all three two-particle cuts gives their eq. (4.12):

\[
\boxed{
M_4^{(1)}(1^+,2^+,3^+,4^+)
=2\frac{[12]^2[34]^2}{\langle12\rangle^2\langle34\rangle^2}
\left(
I_4^{1234}[\mu^8]
+I_4^{3124}[\mu^8]
+I_4^{2314}[\mu^8]
\right).
}
\]

In particular, at four points the all-plus gravity amplitude is **box only** after this D-dimensional reconstruction. No triangle or bubble integral appears in the complete representation. This is a helicity- and multiplicity-specific result, not a generic no-triangle assertion for pure Einstein gravity.

Each box is a scalar box with four massive four-dimensional propagators and numerator `mu^8`.

## Celestial fixed-radius geometry

For a channel of invariant mass `M`, decompose the D-dimensional cut momentum as usual,

\[
\mu=\frac{M}{2\cosh r},
\qquad
\beta=\sqrt{1-\frac{4\mu^2}{M^2}}=\tanh r.
\]

A single `mu^8` box contribution therefore has the universal phase-space radial factor

\[
\boxed{
W_{\rm grav}(r)
=\tanh r\,\operatorname{sech}^8 r,
}
\]

up to the overall `M^8/2^8` and the external graviton-helicity coefficient.

Because a box cut leaves two propagators uncut, its angular dependence is exactly the already-derived massive `S^2` two-propagator master kernel

\[
\mathcal J(r)=
\int_{S^2}\frac{d\Omega}{D_1D_2},
\]

with the closed Feynman-parameter/logarithmic expression recorded in `MASSIVE_S2_TWO_PROPAGATOR_MASTER_KERNEL.md`.

Hence each all-plus gravity box cut has the celestial structure

\[
\boxed{
\text{cut}_{\rm grav}^{(\mu^8)}
\propto
\tanh r\,\sech^8r\;\mathcal J(r)
}
\]

with the appropriate squared external helicity phase.

## Exact normalized radial law

The general `mu^(2k)` shell family at `k=4` gives

\[
\int_0^\infty \tanh r\,\sech^8r\,dr=\frac18.
\]

Therefore

\[
\boxed{
\rho_{\rm grav}(r)=8\tanh r\,\sech^8r
}
\]

is normalized, with

\[
F(R)=1-\sech^8R.
\]

Equivalently,

\[
\boxed{
U=\sech^8r=\left(\frac{2\mu}{M}\right)^8\sim{\rm Uniform}(0,1).
}
\]

The induced transverse-mass density is

\[
\boxed{
\rho_\mu(\mu)=8\left(\frac2M\right)^8\mu^7
=\frac{2048\mu^7}{M^8},
\qquad 0\le\mu\le\frac M2.
}
\]

and

\[
\boxed{
\mathbb E[\mu^q]
=\frac8{q+8}\left(\frac M2\right)^q,
\qquad q>-8.
}
\]

In particular,

\[
\boxed{
\mathbb E\!\left[\left(\frac{2\mu}{M}\right)^2\right]=\frac45.
}
\]

## Unique shell maximum

For `k=4`, the universal maximum condition gives

\[
\boxed{
\tanh r_* = \frac13,
\qquad
\sech^2r_* = \frac89,
\qquad
\mu_* = \frac{\sqrt2}{3}M.
}
\]

The unnormalized maximum is

\[
\boxed{
W_{\rm grav}^{\max}
=\frac13\left(\frac89\right)^4.
}
\]

Near threshold `W_grav(r)=r+O(r^3)`; at the massless boundary

\[
W_{\rm grav}(r)=256e^{-8r}+O(e^{-10r}).
\]

Thus the gravity rational box is much more strongly localized away from the massless boundary than the YM `mu^4` box (`e^{-4r}`) or `mu^2` bubble (`e^{-2r}`).

## Dimension-shifted finite value

Appendix D of Bern et al. gives the required gravity box directly as a 12-dimensional scalar box,

\[
I_4[\mu^8]
=-\epsilon(1-\epsilon)(2-\epsilon)(3-\epsilon)(4\pi)^4
I_4^{12-2\epsilon}.
\]

Its four-dimensional limit is a finite quadratic polynomial in the external invariants. For the massless box it is obtained from their eq. (D.9) by setting the two external masses to zero. This is the mechanism by which the evanescent `mu^8` cut leaves a finite rational four-graviton amplitude.

## Boundary

This closes the D-dimensional topology and celestial radial geometry for the four-point all-plus gravity amplitude. It does **not** imply box-only structure for generic gravity amplitudes or for higher multiplicity: Bern et al. explicitly require pentagons at five and six points, and generic gravity power counting does not exclude triangles/bubbles. It also does not identify the internal hyperbolic radius with the external celestial Mellin scale.
