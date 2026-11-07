# Massive-vector cut: Sym^2 little-group state-sum reduction

Codex/GPT discovery track, 2026-08-25.

The integer-dimensional reconstruction route reduces the nonzero-mu `D_s=4` gauge cut to a four-dimensional massive-vector cut minus one real adjoint-scalar cut. The remaining vector sewing is finite-dimensional because a 4D massive spin-1 state transforms as `Sym^2` of the SU(2) massive little-group doublet.

For the parity-paired all-minus/all-plus massive-vector trees, isolate the two little-group bracket matrices

\[
A_{IJ}=\langle \ell_1^I\,\ell_2^J\rangle,
\qquad
B^{JI}=[\ell_{2,J}\,\ell_{1,I}],
\]

where the reversed ordering in the square bracket is the convention in which the contracted trace is the ordinary Lorentz bispinor trace. Let `C=AB`.

The spin-one tree numerator is the symmetric square of the corresponding spin-one-half matrix. Therefore the sum over the three physical vector states is the trace of the induced endomorphism on `Sym^2(C^2)`:

\[
\boxed{
\operatorname{Tr}_{\mathrm{Sym}^2}(C)
=\frac12\big[(\operatorname{tr}C)^2+\operatorname{tr}(C^2)\big]
=(\operatorname{tr}C)^2-\det C.
}
\]

The second equality is the 2 by 2 Cayley-Hamilton identity.

## Covariant Lorentz reduction

Use the standard massive-spinor completeness relation

\[
\ell_{\alpha\dot\alpha}
=\lambda_\alpha^I\widetilde\lambda_{\dot\alpha I}.
\]

Then the contracted product of the angle and square matrices is

\[
\operatorname{tr}C
=\langle \ell_1^I\,\ell_2^J\rangle
 [\ell_{2,J}\,\ell_{1,I}].
\]

Contracting the little-group indices first reconstructs the two bispinors. The standard two-spinor trace identity gives

\[
\boxed{
\operatorname{tr}C=2\,\ell_1\!\cdot\!\ell_2
}
\]

up to the simultaneous crossing sign attached to both massive legs. Since only `(tr C)^2` enters the spin-one trace, that harmless convention sign drops out.

The determinant is convention invariant. In matrix notation,

\[
A=\lambda_1^{T}\varepsilon\lambda_2,
\]

so

\[
\det A=(\det\lambda_1)(\det\lambda_2)\det\varepsilon.
\]

For an on-shell massive momentum, `det lambda_i` is the mass times a little-group phase/sign fixed by the normalization convention; the conjugate square-bracket matrix carries the inverse phase. Consequently

\[
\boxed{
\det C=(\det A)(\det B)=\mu^4.
}
\]

Hence the complete spin-one little-group contraction is

\[
\boxed{
\operatorname{Tr}_{\mathrm{Sym}^2}(C)
=(2\ell_1\!\cdot\!\ell_2)^2-\mu^4.
}
\]

For an `s`-channel pair with

\[
\ell_1^2=\ell_2^2=\mu^2,
\qquad
(\ell_1+\ell_2)^2=s,
\]

we have `2 ell_1 dot ell_2=s-2 mu^2`, and therefore

\[
\boxed{
\operatorname{Tr}_{\mathrm{Sym}^2}(C)
=(s-2\mu^2)^2-\mu^4
=s^2-4s\mu^2+3\mu^4.
}
\]

At threshold `s=4 mu^2`, this becomes `3 mu^4`, exactly matching the three physical massive-vector polarizations. The companion `audit_massive_vector_sym2.py` checks the same formula in an explicit center-of-mass massive-spinor frame.

## Dimensional-reconstruction consequence

The state-counting identity remains

\[
C^{(D_s=4)}(\mu)=C^{(V_m)}(\mu)-C^{(S)}(\mu).
\]

The vector polarization numerator is therefore no longer an unspecified state sum: its universal little-group factor is the polynomial above. What remains before declaring the complete `D_s=4` MHV cut closed is to insert the exact massive-vector tree denominators/external helicity factors, align their normalization with the real-adjoint-scalar convention, and perform the subtraction. That normalization check is essential because the scalar term may be packaged as one real or one complex species in amplitude references.
