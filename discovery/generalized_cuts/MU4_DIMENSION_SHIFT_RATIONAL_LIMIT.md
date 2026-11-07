# mu^4 scalar box: dimension-shift rational limit

Codex/GPT discovery track, 2026-08-25.

## Normalization

Use the standard dimensionally regulated scalar-integral convention

\[
I_4^D:=\int\frac{d^D L}{i\pi^{D/2}}
\frac{1}{D_1D_2D_3D_4},
\]

with `mu^2 = -L_perp^2` in the split `D=4-2 epsilon` loop momentum. Overall factors such as
`(4 pi)^epsilon`, `e^{-gamma_E epsilon}`, coupling constants, and the amplitude-specific
helicity numerator are not included here.

## Exact dimension shift

The standard transverse-moment identity is

\[
I_n^{4-2\epsilon}[\mu^{2r}]
=\frac{\Gamma(r-\epsilon)}{\Gamma(-\epsilon)}
I_n^{4+2r-2\epsilon}.
\]

For `r=2`, Gamma recursion gives

\[
\frac{\Gamma(2-\epsilon)}{\Gamma(-\epsilon)}
=(-\epsilon)(1-\epsilon),
\]

hence

\[
\boxed{
I_4^{4-2\epsilon}[\mu^4]
=-\epsilon(1-\epsilon) I_4^{8-2\epsilon}.
}
\]

## Universal UV residue of the eight-dimensional box

Feynman parameterization gives

\[
I_4^D
=\Gamma\!\left(4-\frac D2\right)
\int_{x_i\ge0}\!d^4x\,\delta\!\left(1-\sum_i x_i\right)
[\mathcal F(x)-i0]^{D/2-4}.
\]

At `D=8-2 epsilon`,

\[
I_4^{8-2\epsilon}
=\Gamma(\epsilon)
\int_{\Delta_3} [\mathcal F(x)-i0]^{-\epsilon}\,d^3x.
\]

For generic nonexceptional external kinematics the parameter integral tends to the volume of
the standard 3-simplex,

\[
\operatorname{Vol}(\Delta_3)=\frac1{3!}=\frac16.
\]

Since

\[
\Gamma(\epsilon)=\frac1\epsilon+O(1),
\]

we obtain the universal pole

\[
\boxed{
I_4^{8-2\epsilon}=\frac{1}{6\epsilon}+O(1).
}
\]

Therefore

\[
I_4^{4-2\epsilon}[\mu^4]
=-\epsilon(1-\epsilon)
\left(\frac{1}{6\epsilon}+O(1)\right),
\]

and the regulator limit closes to

\[
\boxed{
\lim_{\epsilon\to0}
I_4^{4-2\epsilon}[\mu^4]
=-\frac16.
}
\]

This explains precisely how the `mu^4` numerator can vanish pointwise on a strict 4D cut yet
leave a finite rational term after D-dimensional integration: the explicit `epsilon` from the
dimension shift multiplies the UV `1/epsilon` pole of the dimension-raised box.

## Relation to the celestial radial shell

The fixed-radius celestial cut isolated earlier has `mu^4 = M^4 sech^4(r)/16` and therefore
vanishes at the massless boundary `r -> infinity`. The finite `-1/6` result is not a contradiction:
it is produced only after integrating over the full D-dimensional transverse sector and taking
the dimension-shifted UV residue. Thus the radial shell description and the rational term are two
representations of the same evanescent information at different stages of the calculation.

## Boundary

The constant `-1/6` is the scalar integral insertion in the normalization stated above. It is not
yet the complete pure-Yang--Mills rational amplitude coefficient. The external helicity phase,
color/coupling normalization, the `D_s=4, mu!=0` gluon baseline, scalar-species convention, and
triangle/bubble sectors must still be assembled consistently.
