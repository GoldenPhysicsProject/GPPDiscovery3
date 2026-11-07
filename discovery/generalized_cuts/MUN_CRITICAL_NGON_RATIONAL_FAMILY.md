# Critical mu-power n-gon rational family

Codex/GPT discovery track, 2026-08-25.

Use the same scalar-integral normalization

\[
I_n^D:=\int\frac{d^D L}{i\pi^{D/2}}\prod_{j=1}^n\frac1{D_j},
\]

and the standard transverse-moment dimension shift

\[
I_n^{4-2\epsilon}[\mu^{2r}]
=\frac{\Gamma(r-\epsilon)}{\Gamma(-\epsilon)}
I_n^{4+2r-2\epsilon}.
\]

For integer `r>=1`, Gamma recursion gives

\[
\frac{\Gamma(r-\epsilon)}{\Gamma(-\epsilon)}
=(-\epsilon)(1-\epsilon)(2-\epsilon)\cdots(r-1-\epsilon)
=-(r-1)!\,\epsilon+O(\epsilon^2).
\]

Now choose the critical power

\[
\boxed{r=n-2.}
\]

Then the raised dimension is

\[
D'=4+2(n-2)-2\epsilon=2n-2\epsilon.
\]

The Feynman-parameter prefactor of an n-gon is

\[
\Gamma\!\left(n-\frac{D'}2\right)=\Gamma(\epsilon)
=\frac1\epsilon+O(1).
\]

At `epsilon=0` the parameter exponent is zero, so the remaining parameter integral tends, for generic nonexceptional kinematics, to the standard `(n-1)`-simplex volume

\[
\operatorname{Vol}(\Delta_{n-1})=\frac1{(n-1)!}.
\]

Hence

\[
I_n^{2n-2\epsilon}
=\frac1{(n-1)!\,\epsilon}+O(1).
\]

Multiplying by the dimension-shift factor gives the universal finite limit

\[
\boxed{
\lim_{\epsilon\to0}
I_n^{4-2\epsilon}[\mu^{2(n-2)}]
=-\frac{(n-3)!}{(n-1)!}
=-\frac1{(n-1)(n-2)}
}
\]

for every integer `n>=3`, in the stated normalization.

Examples:

\[
n=3:\quad I_3[\mu^2]\to-\frac12,
\]

\[
n=4:\quad I_4[\mu^4]\to-\frac16,
\]

\[
n=5:\quad I_5[\mu^6]\to-\frac1{12},
\]

\[
n=6:\quad I_6[\mu^8]\to-\frac1{20}.
\]

This gives a universal hierarchy of evanescent rational residues. The box result `-1/6` is the `n=4` member, not an isolated coincidence.

The mechanism is always the same:

1. the critical `mu` numerator vanishes in a strict four-dimensional projection;
2. dimensional shifting supplies one explicit factor of `epsilon`;
3. the raised `2n`-dimensional n-gon supplies the matching logarithmic UV `1/epsilon` pole;
4. the simplex volume fixes the rational remainder.

Boundary: this is a scalar-integral insertion identity in the stated normalization. It does not imply that a given Yang--Mills or gravity amplitude contains the critical numerator with unit coefficient, nor that lower topologies or state-sum terms vanish. Those coefficients must come from honest D-dimensional generalized cuts.
