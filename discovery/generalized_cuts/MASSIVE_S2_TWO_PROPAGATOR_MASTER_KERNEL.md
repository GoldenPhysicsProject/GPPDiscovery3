# Massive fixed-radius S2 cut: exact two-propagator master kernel

Codex/GPT discovery track, 2026-08-24.

## Scope

The D-dimensional MHV massive-scalar sector has the form

\[
C_s^{\rm scalar}=\mu^4\,\Xi(1,2,3,4)\,\frac1{D_1^{(\mu)}D_2^{(\mu)}}.
\]

At fixed channel mass `M` and transverse mass `mu`, the two-particle cut is the angular sphere `S^2` at one fixed radius in `H^3`. This note performs the universal angular reduction of the two uncut propagators. The external helicity factor `Xi` is left untouched.

## Universal affine form on the cut sphere

In the channel center-of-mass frame, write the cut momentum as

\[
\ell=(E,p\,\hat n),\qquad \hat n\in S^2,
\]

with

\[
E=\frac M2,\qquad p=\frac12\sqrt{M^2-4\mu^2}.
\]

Any uncut massive propagator

\[
D^{(\mu)}(K)=(\ell+K)^2-\mu^2
\]

is affine in `hat n`:

\[
D^{(\mu)}(K)=A_K+B_K\cdot\hat n,
\]

where, for metric `(+---)`,

\[
A_K=K^2+2E K^0,
\qquad
B_K=-2p\,\mathbf K.
\]

Thus the entire angular problem reduces to

\[
\mathcal J(A,B;C,D)
:=\int_{S^2}\frac{d\Omega(\hat n)}{(A+B\cdot\hat n)(C+D\cdot\hat n)}.
\]

The formulas below hold in the nonsingular Euclidean angular domain in which the affine denominators do not vanish on `S^2`; analytic continuation gives the physical `i0` prescription.

## Exact Feynman-parameter reduction

Use

\[
\frac1{XY}=\int_0^1\frac{dx}{[xX+(1-x)Y]^2}.
\]

Define

\[
\alpha(x)=xA+(1-x)C,
\qquad
V(x)=xB+(1-x)D.
\]

Rotational invariance gives the elementary sphere integral

\[
\int_{S^2}\frac{d\Omega}{(\alpha+V\cdot\hat n)^2}
=\frac{4\pi}{\alpha^2-|V|^2}.
\]

Therefore

\[
\boxed{
\mathcal J(A,B;C,D)
=4\pi\int_0^1\frac{dx}{Q(x)},
}
\]

with the quadratic

\[
Q(x)=\alpha(x)^2-|V(x)|^2=q_2x^2+q_1x+q_0,
\]

and exact coefficients

\[
\boxed{q_0=C^2-|D|^2,}
\]

\[
\boxed{q_1=2\,[C(A-C)-D\cdot(B-D)],}
\]

\[
\boxed{q_2=(A-C)^2-|B-D|^2.}
\]

This is already a complete angular reduction: no `S^2` integration remains.

## Closed logarithmic form

Let

\[
\Delta=q_1^2-4q_2q_0.
\]

For `q2 != 0` and a branch of `sqrt(Delta)` chosen consistently with the `i0` prescription,

\[
\int\frac{dx}{q_2x^2+q_1x+q_0}
=\frac1{\sqrt\Delta}
\log\frac{2q_2x+q_1-\sqrt\Delta}{2q_2x+q_1+\sqrt\Delta}.
\]

Hence

\[
\boxed{
\mathcal J
=\frac{4\pi}{\sqrt\Delta}
\log\left[
\frac{(2q_2+q_1-\sqrt\Delta)(q_1+\sqrt\Delta)}
{(2q_2+q_1+\sqrt\Delta)(q_1-\sqrt\Delta)}
\right].
}
\]

For `Delta<0` this is equivalently an arctangent form after combining conjugate logarithms. Degenerate `q2=0` or `Delta=0` cases are obtained by the elementary linear/repeated-root limits.

## Massive-scalar MHV sector

The fixed-radius scalar-sector cut therefore becomes

\[
\boxed{
\int d\Pi_2\,C_s^{\rm scalar}
=\frac{\beta}{32\pi^2}\,
\mu^4\,\Xi\,\mathcal J,
\qquad
\beta=\sqrt{1-\frac{4\mu^2}{M^2}},
}
\]

when `dPi2 = beta dOmega/(32 pi^2)` is used. Since

\[
\mu^4=\frac{M^4}{16\cosh^4 r},
\qquad
\beta=\tanh r,
\]

this is equivalently

\[
\boxed{
\int d\Pi_2\,C_s^{\rm scalar}
=\frac{M^4}{512\pi^2}
\frac{\tanh r}{\cosh^4 r}\,\Xi\,\mathcal J(r).
}
\]

Thus the D-dimensional rational-box numerator contributes the exact radial prefactor

\[
\boxed{\tanh r\,\operatorname{sech}^4 r}
\]

multiplying a logarithmic angular master kernel. This is stronger than the earlier observation of a bare `sech^4 r`: inclusion of the two-body measure supplies the additional `tanh r` threshold factor.

## Threshold and massless radial behavior

At threshold `r -> 0`,

\[
\tanh r\,\operatorname{sech}^4r=r+O(r^3),
\]

so this sector is phase-space suppressed linearly, modulo nonsingular behavior of the propagator kernel.

At the massless boundary `r -> infinity`,

\[
\tanh r\,\operatorname{sech}^4r\sim16e^{-4r},
\]

which matches the fact that the isolated `mu^4` numerator vanishes pointwise in a strict four-dimensional cut. Its finite rational contribution can therefore arise only after the full D-dimensional integral/dimension-shift mechanism is retained; one must not infer vanishing of the integrated rational term from this fixed-`mu` boundary behavior.

## Boundary

This closes the universal `S^2` angular integral for the two uncut massive propagators but does not yet fix the exact cyclic spinor phase `Xi`, assemble the full D-dimensional gluon state sum, perform triangle/bubble subtraction, or execute the dimension-shifted loop integration. Those are the remaining requirements for an honest pure-Yang--Mills rational amplitude. Gravity remains downstream of that state-sum closure.
