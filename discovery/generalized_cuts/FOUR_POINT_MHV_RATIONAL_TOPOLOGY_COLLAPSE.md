# Four-point adjacent-MHV rational sector collapses to mu^4 box plus mu^2 bubble

Codex/GPT discovery track, 2026-08-25.

## External input checked against Badger 0806.4600

Badger's D-dimensional generalized-unitarity basis writes the rational part in terms of

\[
C_4^{[4]} I_4[\mu^4],\qquad
C_3^{[2]} I_3[\mu^2],\qquad
C_2^{[2]} I_2[\mu^2].
\]

For the four-point adjacent-MHV helicity ordering

\[
1^-\,2^-\,3^+\,4^+,
\]

his explicit four-point coefficient list gives the box coefficient

\[
\boxed{
C_4^{[4]}(1^-,2^-,3^+,4^+)
=2i\,\frac{\langle12\rangle[43]}
          {\langle34\rangle[21]}
=2i\,\frac{\langle12\rangle[34]}
          {\langle34\rangle[12]}.
}
\]

This agrees with the convention-fixed scalar-cut phase already obtained here,

\[
\Xi=-\frac{\langle12\rangle[34]}{[12]\langle34\rangle},
\]

through

\[
\boxed{C_4^{[4]}=-2i\,\Xi.}
\]

The same explicit coefficient list gives triangle `mu^2` coefficients for the one-minus and alternating-MHV helicity sectors and then states that all other triangle coefficients vanish. Therefore, for the adjacent-MHV ordering above,

\[
\boxed{C_3^{[2]}(1^-,2^-,3^+,4^+)=0}
\]

for every triangle channel.

The bubble list contains a nonzero adjacent-MHV coefficient (Badger eq. 5.33, in the `23` channel); all other bubble coefficients not related by the cyclic/complement identification vanish.

## Consequence

Combining this with the supersymmetric rational reduction gives the complete topology content required for the four-point pure-Yang--Mills rational remainder in the FDH organization:

\[
\boxed{
R_4^{\rm YM}(--++)
\;\Longleftarrow\;
C_4^{[4]} I_4[\mu^4]
\; +\;
C_2^{[2]} I_2[\mu^2],
}
\]

with **no `mu^2` triangle contribution**.

The box side is already under control in the present workbench:

- the massive-scalar helicity factor is fixed;
- the fixed-radius `S^2` two-propagator angular kernel is closed;
- the dimension-shift limit is `I_4[mu^4] -> -1/6` in the stated scalar-integral normalization.

Hence the sole unresolved lower-topology generalized-cut object for this four-point MHV rational sector is the scalar `mu^2` bubble coefficient.

## Bubble dimension-shift normalization

Badger's general rational formula weights each independent bubble coefficient by the finite rational contribution of `I_2[mu^2]`; in the massless two-particle channel this is proportional to the channel invariant. The next calculation should therefore derive the bubble coefficient directly from the massive scalar double cut and keep the channel invariant explicit, rather than importing the final amplitude.

## Boundary

This statement is specific to the four-point adjacent-MHV scalar rational sector and the FDH/supersymmetric decomposition. It does not say triangle coefficients vanish at higher multiplicity or for other helicity orderings, and it does not yet include the finite FDH-to-HV scheme conversion.

Source checked directly: S. D. Badger, JHEP 01 (2009) 049, arXiv:0806.4600, especially eqs. (5.1)--(5.3), (5.10), the triangle list (5.17)--(5.24), and the bubble list including eq. (5.33).
