# D-dimensional MHV cut: explicit massive-scalar mu^4 sector

Codex/GPT discovery track, 2026-08-24.

## Scope

A D-dimensional massless loop momentum decomposes as

\[
L=\ell+L_\perp,\qquad \mu^2=-L_\perp^2,
\]

so the four-dimensional component obeys

\[
\ell^2=\mu^2.
\]

Thus the D-dimensional cut can be evaluated with four-dimensional **massive** cut legs. The companion H3 note identifies their fixed-radius phase-space geometry. Here we isolate the first numerator information that disappears if one sets `mu=0` too early.

## Massive scalar tree input

Badger, *Direct Extraction Of One Loop Rational Terms* (arXiv:0806.4600), gives the color-ordered tree amplitudes

\[
A_4(1_s,2^+,3^+,4_s)
=i\,\frac{\mu^2[23]}
{\langle23\rangle\langle2|1|2]},
\]

and

\[
A_4(1_s,2^+,3^-,4_s)
=i\,\frac{\langle3|1|2]^2}
{s_{23}\langle2|1|2]}.
\]

Here legs `1_s,4_s` are massive scalars of mass `mu`. Parity gives

\[
A_4(1_s,2^-,3^-,4_s)
=i\,\frac{\mu^2\langle23\rangle}
{[23]\,[2|1|2\rangle},
\]

up to the same all-outgoing spinor-phase conventions used for the conjugate amplitude.

The crucial invariant fact is independent of those phases:

\[
A_4(s,--,s)\propto\mu^2,
\qquad
A_4(s,++,s)\propto\mu^2.
\]

## The scalar sector of the s-channel MHV cut

For external helicities

\[
1^-\,2^-\,3^+\,4^+,
\]

use the same all-outgoing routing as the four-dimensional gluon-cut note,

\[
p_1+p_2+\ell_1+\ell_2=0.
\]

The adjoint-scalar internal-state contribution is the product

\[
C_s^{\rm scalar}
=
A_4(\ell_{1,s},1^-,2^-,\ell_{2,s})
A_4((-\ell_2)_s,3^+,4^+,(-\ell_1)_s),
\]

with the precise cyclic ordering adjusted to the chosen color routing.

Each tree supplies one factor `mu^2`, hence

\[
\boxed{C_s^{\rm scalar}\propto\mu^4.}
\]

The spinor-sandwich denominators are just uncut massive propagators. For a massive cut leg `ell^2=mu^2` and a massless external momentum `p`,

\[
(\ell+p)^2-\mu^2=2\ell\cdot p
=\langle p|\ell|p],
\]

up to the standard bracket orientation convention. Therefore the scalar-sector cut has the structural form

\[
\boxed{
C_s^{\rm scalar}
=
\mu^4\,\Xi(1,2,3,4)\,
\frac{1}{D_1^{(\mu)}D_2^{(\mu)}}
}
\]

where `Xi` is an external spinor phase/rational helicity factor and

\[
D_i^{(\mu)}=(\ell+K_i)^2-\mu^2
\]

are the two uncut massive scalar propagators.

This is already enough to establish the key physics point: **the D-dimensional MHV cut contains a box-sector numerator proportional to `mu^4` which identically vanishes on a strictly four-dimensional cut.**

## Why this is exactly the rational-term sector one expects

Badger's D-dimensional generalized-unitarity analysis proves the power-counting bound that box coefficients in a renormalizable gauge theory can contain terms through `mu^4`, whereas triangle and bubble coefficients require at most `mu^2`. The rational box coefficient is extracted precisely from the `mu^4` boundary term of the product of massive tree amplitudes.

Thus the `mu^4` product above is not an arbitrary regulator correction: it lies in the exact highest-dimensional box sector responsible for rational information missed by ordinary four-dimensional unitarity.

## Relation to the H3 geometry

At fixed `mu` and channel mass `M`, each cut leg lies on the H3 shell

\[
\cosh r=\frac{M}{2\mu},
\qquad
\tanh r=\sqrt{1-\frac{4\mu^2}{M^2}}.
\]

Therefore the D-dimensional numerator and celestial geometry now meet concretely:

- the angular cut is the `S^2` at fixed hyperbolic radius `r`;
- the massive scalar state sector contributes `mu^4` times two uncut massive propagators;
- `mu=M/(2\cosh r)`, so the missing numerator can equivalently be written as

\[
\boxed{
\mu^4=\frac{M^4}{16\cosh^4 r}.
}
\]

This converts the rational-term numerator into an explicit radial weight on the massive celestial hyperboloid.

## Immediate analytic target

The next calculation is to insert the explicit scalar-tree product, including its external helicity phase, into the fixed-radius H3/S2 cut integral and perform the angular reduction. That will determine the massive analogue of the scalar-box cut kernel with the new radial weight `cosh^{-4} r`.

Only after that reduction is it legitimate to compare the resulting `mu^4` contribution with the already-controlled scalar regulator dispersion theorem.

## Boundary

This note isolates the massive-scalar sector used in D-dimensional unitarity. It does not yet assemble the complete scheme-dependent D-dimensional gluon state sum, triangle/bubble subtraction terms, or the final pure-Yang--Mills rational amplitude. It also does not identify the earlier infrared regulator `m` with this transverse `mu^2`.

Reference: S. D. Badger, arXiv:0806.4600, especially the massive scalar tree amplitudes (eqs. 56--57) and the `mu^4` box-coefficient extraction (eqs. 32--33 and Appendix A).
