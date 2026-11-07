# First-principles simplex derivation of the massless mu^8 gravity box residue

Codex/GPT discovery track, 2026-08-25.

The Appendix-D massless gravity-box polynomial

\[
I_4[\mu^8]\propto
-\frac{2s^2+2t^2+st}{840}
\]

can be derived directly from the dimension shift and elementary moments of the Feynman-parameter simplex. This removes the need to treat the `1/840` as an imported lookup constant.

## 1. Dimension shift

For the scalar box in `D=4-2epsilon`, insertion of `mu^8` means `r=4`:

\[
I_4^{4-2\epsilon}[\mu^8]
=
\frac{\Gamma(4-\epsilon)}{\Gamma(-\epsilon)}
I_4^{12-2\epsilon}.
\]

Gamma recursion gives

\[
\frac{\Gamma(4-\epsilon)}{\Gamma(-\epsilon)}
=(-\epsilon)(1-\epsilon)(2-\epsilon)(3-\epsilon)
=-6\epsilon+O(\epsilon^2).
\]

The raised-dimensional box has

\[
4-\frac{12-2\epsilon}{2}=-2+\epsilon,
\]

so its Feynman-parameter prefactor is

\[
\Gamma(-2+\epsilon)
=\frac1{2\epsilon}+O(1).
\]

Therefore the product of the evanescent dimension-shift factor and the ultraviolet pole tends to

\[
\boxed{-3}.
\]

## 2. Squared Feynman polynomial

For a massless box with adjacent invariants `s,t`, the Symanzik/Feynman polynomial is, up to the conventional overall sign irrelevant here because it is squared,

\[
F=s\,a_1a_3+t\,a_2a_4.
\]

At `D=12`, the parameter exponent is

\[
D/2-n=6-4=2,
\]

so the residue requires

\[
\int_{\Delta_3}F^2.
\]

Expanding,

\[
F^2
=s^2a_1^2a_3^2
+t^2a_2^2a_4^2
+2st\,a_1a_2a_3a_4.
\]

## 3. Dirichlet simplex moments

For the standard three-simplex

\[
\Delta_3=\{a_i\ge0:\ a_1+a_2+a_3+a_4=1\},
\]

the elementary Dirichlet moment formula is

\[
\int_{\Delta_3}\prod_{i=1}^4 a_i^{m_i}
=\frac{\prod_i m_i!}{(3+\sum_i m_i)!}.
\]

Hence

\[
\int_{\Delta_3}a_1^2a_3^2
=\frac{2!\,2!}{7!}
=\frac1{1260},
\]

\[
\int_{\Delta_3}a_2^2a_4^2=\frac1{1260},
\]

and

\[
\int_{\Delta_3}a_1a_2a_3a_4
=\frac1{7!}
=\frac1{5040}.
\]

Therefore

\[
\begin{aligned}
\int_{\Delta_3}F^2
&=\frac{s^2+t^2}{1260}+\frac{2st}{5040}\\
&=\boxed{\frac{2s^2+2t^2+st}{2520}}.
\end{aligned}
\]

## 4. Finite rational residue

Multiplying by the dimension-shift/pole factor `-3` gives

\[
\boxed{
I_4[\mu^8]
\longrightarrow
-\frac{2s^2+2t^2+st}{840}
}
\]

apart from the loop-measure convention factor `i/(4pi)^2` carried separately by Bern et al.

This is exactly the massless specialization of their Appendix-D eq. (D.9).

## Structural interpretation

The earlier critical family `mu^{2(n-2)}` corresponds to parameter exponent zero and therefore produces only the simplex volume. The gravity box is **supercritical**: for `n=4`, `r=4`, the excess

\[
m=r-(n-2)=2
\]

forces the second moment `F^2` of the simplex geometry to survive. That is why the Yang--Mills `mu^4` box residue is kinematics-independent while the gravity `mu^8` box residue is a quadratic polynomial in Mandelstam invariants.

This suggests the general supercritical pattern: an insertion `mu^{2r}` with `r=n-2+m` probes the `m`th Feynman-polynomial moment of the `(n-1)`-simplex.
