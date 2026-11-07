# Integrated four-graviton all-plus closure from the three mu^8 boxes

Codex/GPT discovery track, 2026-08-25.

This completes the four-dimensional limit of the box representation in `FOUR_GRAVITON_ALL_PLUS_MU8_BOX_CLOSURE.md` directly from Appendix D of Bern--Dixon--Perelstein--Rozowsky.

## Massless limit of the gravity box

Their Appendix-D formula (D.9) gives, after setting both external masses to zero,

\[
I_4[\mu^8](x,y)
=-\frac{i}{(4\pi)^2}\frac{1}{840}
\left(2x^2+2y^2+xy\right)+O(\epsilon),
\]

where `x,y` are the two adjacent Mandelstam invariants of the massless box.

Define

\[
P(x,y):=2x^2+2y^2+xy.
\]

For the three box orderings in eq. (4.12),

\[
I_4^{1234}[\mu^8]\leftrightarrow P(s,t),
\]

\[
I_4^{3124}[\mu^8]\leftrightarrow P(u,s),
\]

\[
I_4^{2314}[\mu^8]\leftrightarrow P(t,u),
\]

with

\[
s=s_{12},\qquad t=s_{23},\qquad u=s_{13},\qquad s+t+u=0.
\]

## Exact polynomial collapse

The sum is

\[
P(s,t)+P(u,s)+P(t,u)
=4(s^2+t^2+u^2)+(st+tu+us).
\]

Since

\[
(s+t+u)^2=0
\quad\Longrightarrow\quad
st+tu+us=-\frac12(s^2+t^2+u^2),
\]

we obtain

\[
\boxed{
P(s,t)+P(u,s)+P(t,u)
=\frac72(s^2+t^2+u^2).
}
\]

Therefore

\[
\sum_{\rm 3\ boxes} I_4[\mu^8]
=-\frac{i}{(4\pi)^2}
\frac{s^2+t^2+u^2}{240}.
\]

The complete amplitude representation carries an additional overall factor of two,

\[
M_4^{(1)}(++++)
=2\frac{[12]^2[34]^2}
        {\langle12\rangle^2\langle34\rangle^2}
\sum_{\rm 3\ boxes}I_4[\mu^8].
\]

Hence

\[
\boxed{
M_4^{(1)}(++++)
=-\frac{i}{(4\pi)^2}
\frac{[12]^2[34]^2}
     {\langle12\rangle^2\langle34\rangle^2}
\frac{s^2+t^2+u^2}{120}.
}
\]

Using the standard four-point spinor identity

\[
\frac{[12]^2[34]^2}
     {\langle12\rangle^2\langle34\rangle^2}
=
\left(
\frac{st}
{\langle12\rangle\langle23\rangle\langle34\rangle\langle41\rangle}
\right)^2,
\]

this is exactly Bern et al. eq. (4.14):

\[
\boxed{
M_4^{(1)}(1^+,2^+,3^+,4^+)
=-\frac{i}{(4\pi)^2}
\left(
\frac{st}
{\langle12\rangle\langle23\rangle\langle34\rangle\langle41\rangle}
\right)^2
\frac{s^2+t^2+u^2}{120}.
}
\]

## Structural consequence

The entire finite all-plus four-graviton rational amplitude is therefore the integrated image of three `mu^8` hyperbolic box shells. The nontrivial symmetric polynomial `s^2+t^2+u^2` arises from summing the three channel-dependent quadratic Appendix-D residues, not from assigning one universal constant to each gravity box.

This is an important difference from the Yang--Mills `mu^4` box insertion, whose critical raised-dimensional box residue is kinematics-independent.

Primary source: arXiv:hep-th/9811140, eqs. (4.12)--(4.14), definitions below (D.4), and eq. (D.9).
