# Massive-scalar MHV cut: closure of the external helicity factor Xi

Codex/GPT discovery track, 2026-08-25.

Use exactly the cyclic routing already adopted for the scalar contribution to the
`1^- 2^- 3^+ 4^+` s-channel cut:

\[
C_s^{\rm scalar}
=A_4(\ell_{1,s},1^-,2^-,\ell_{2,s})
 A_4((-\ell_2)_s,3^+,4^+,(-\ell_1)_s).
\]

From the massive-scalar tree formulas recorded in `D_DIMENSIONAL_MHV_MASSIVE_SCALAR_SECTOR.md`,

\[
A_4(1_s,2^-,3^-,4_s)
=i\mu^2\frac{\langle23\rangle}{[23]\,[2|1|2\rangle},
\]

and

\[
A_4(1_s,2^+,3^+,4_s)
=i\mu^2\frac{[23]}{\langle23\rangle\,\langle2|1|2]}.
\]

Therefore, defining the two uncut denominator sandwiches by

\[
D_L=[1|\ell_1|1\rangle,
\qquad
D_R=\langle3|(-\ell_2)|3],
\]

the two trees are

\[
A_L=i\mu^2\frac{\langle12\rangle}{[12]D_L},
\qquad
A_R=i\mu^2\frac{[34]}{\langle34\rangle D_R}.
\]

Their product is thus

\[
\boxed{
C_s^{\rm scalar}
=-\mu^4
\frac{\langle12\rangle[34]}{[12]\langle34\rangle}
\frac1{D_LD_R}.
}
\]

Hence, in this explicit sandwich convention,

\[
\boxed{
\Xi(1,2,3,4)
=-\frac{\langle12\rangle[34]}{[12]\langle34\rangle}.
}
\]

For real Lorentzian external momenta, square and angle brackets are conjugate up to the standard crossing phases, so this ratio has unit modulus whenever the displayed brackets are nonzero:

\[
|\Xi|=1.
\]

Thus the magnitude of the highest-dimensional scalar-box numerator is exactly `mu^4`; all external helicity information in this sector is carried by a pure little-group phase.

When `D_L,D_R` are rewritten as quadratic propagators `(ell+K)^2-mu^2`, an overall sign can move between `Xi` and the definitions of the `D_i` depending on the all-outgoing routing. The sandwich formula above is therefore the convention-safe closure; any quadratic-propagator version must state its routing explicitly.

This closes the previously abstract `Xi` in the isolated massive-scalar box sector. It still does not assemble the full `D_s` gluon state sum or triangle/bubble subtraction sectors.
