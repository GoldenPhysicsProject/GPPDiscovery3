# Four-point adjacent-MHV bubble: Lean promotion boundary

Codex/GPT discovery track, 2026-08-26.

The direct Badger/Forde scalar-loop audit now closes the `s23` adjacent-MHV
`mu^2` bubble coefficient.  This note separates what is pure algebra and can be
promoted to Lean immediately from what remains a convention-dependent physics input.

## Lean-ready algebra 1: triangle-pole discriminant

Let

\[
B=t(1-u^2)-u,
\qquad
C=u\mu^2-u t^2+u^2t.
\]

For the quadratic

\[
P(y)=u y^2+B y+C,
\]

the discriminant satisfies the polynomial identity

\[
\boxed{
B^2-4uC
=\bigl[t(1+u^2)-u\bigr]^2-4\mu^2u^2.
}
\]

No cut convention, branch choice, or nonzero denominator is required for this identity.
It is a direct `ring` theorem over `ℝ`.

## Lean-ready algebra 2: moment-mapped subtraction

After the exact Vieta sum and Badger `T_1,T_2,T_3` moment map, the symbolic audit
reduces the surviving subtraction to

\[
\boxed{
\mathcal T(u,\mu^2)
=\frac{i\,(10\mu^2u^2+6\mu^2-u^2)}{3(1+u^2)}.
}
\]

Consequently its `mu^2` coefficient is

\[
\boxed{
[\mu^2]\mathcal T
=\frac{2i(5u^2+3)}{3(1+u^2)}.
}
\]

For real `u`, `1+u^2>0`, so the denominator is never singular.  The corresponding
real coefficient identity can be formalized without complex differentiation by
introducing the affine numerator explicitly or simply proving the difference quotient.

After Badger's explicit `-1/2` subtraction prefactor, one scalar-flow orientation is

\[
C_{\rm one}
=-\frac{i(5u^2+3)}{3(1+u^2)}.
\]

## Lean-ready algebra 3: frame restoration conditional on multiplicity two

The rational frame obeys

\[
s_{23}=1,
\qquad
s_{12}=-\frac{u^2}{1+u^2},
\qquad Q=1.
\]

Hence the exact rational identity

\[
\boxed{
2\left[-\frac{5u^2+3}{3(1+u^2)}\right]
=\frac23\left(2\left[-\frac{u^2}{1+u^2}\right]-3\right)
}
\]

can be promoted directly. Multiplying both sides by `i` gives the frame form of the
published bubble coefficient.

## Physics input that must remain explicit

The equality

\[
C_{\rm complex\ scalar}=2 C_{\rm one}
\]

is **not** a consequence of the rational identities above.  It uses the convention
that the `A^[s]` sector is one complex scalar, together with equality of the two
real-component / scalar-flow contributions for purely gluonic external states.

Therefore a proof-assistant theorem that reconstructs the full scalar-loop bubble from
the one-flow result should either:

1. take `Ccomplex = 2 * Cone` as an explicit hypothesis, or
2. separately formalize a state-space model in which the complex scalar is represented
   as two equal real scalar components.

It must not hide this input inside a definition and then present the full coefficient as
though the multiplicity had been derived from the cut algebra.

## Pure-bubble boundary

Independently, the exact double-cut `Inf_y` boundary is

\[
-\frac{tu^2}{1+u^2}\bigl(tu^2+3t+uy-u\bigr),
\]

which contains no `mu^2`.  Therefore the pure-bubble `mu^2` coefficient is exactly zero.
This is also Lean-ready once the polynomial boundary expression itself is supplied as a
hypothesis or represented by a definition; deriving that expression from spinor trees is
still outside the current Lean model.

## Recommended theorem boundary

The next formal file should prove the discriminant identity, denominator positivity,
closed moment-map coefficient identity, and frame restoration.  It should stop there.
The existing `FourPointMHVRationalClosure.lean` can then consume the reduced coefficient
only after the complex-scalar multiplicity is supplied explicitly.
