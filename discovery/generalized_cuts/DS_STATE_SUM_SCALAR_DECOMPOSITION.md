# D_s-dimensional gluon state sum: dimensional reconstruction baseline plus adjoint scalars

Codex/GPT discovery track, 2026-08-25.

For a massless internal gauge boson whose spin states live in `D_s` dimensions, the physical state count is `D_s-2`.  In dimensional reconstruction, however, one must distinguish the spin dimension `D_s` from the loop-momentum dimension `D`.

Keep the same D-dimensional loop momentum

\[
L=\ell+L_\perp,\qquad \ell^2=\mu^2,
\]

and compare state sums at different values of `D_s`.  Increasing `D_s` by one adds one metric direction orthogonal to the fixed four-dimensional external kinematics.  By gauge invariance that extra polarization couples to external four-dimensional gluons as one real adjoint scalar species.  Consequently the state-sum dependence on `D_s` is linear:

\[
\boxed{
C^{(D_s)}(\mu)
=C^{(D_s=4)}(\mu)
 +(D_s-4)\,C^{(\mathrm{one\ real\ adjoint\ scalar})}(\mu).
}
\]

The coefficient is `D_s-4`, not `(D_s-4)^2`, because the scalar flavour index is transported across the two-particle cut and traced once around the loop.

If the scalar-tree convention packages two real adjoint scalars into one complex scalar, then

\[
C^{(D_s)}(\mu)
=C^{(D_s=4)}(\mu)
 +\frac{D_s-4}{2}\,C^{(\mathrm{one\ complex\ adjoint\ scalar})}(\mu).
\]

## Important correction to the first version

The baseline `C^(D_s=4)(mu)` must **not** be identified with the strict four-dimensional massless-helicity cut when `mu != 0`.  The loop momentum is still D-dimensional, so its four-dimensional component is massive.  Dimensional reconstruction changes the spin-state dimension while holding that loop momentum fixed.

Only after separately taking `mu -> 0` does the baseline reduce to the ordinary strict 4D helicity cut already recorded in `YM_GRAVITY_MHV_TWO_PARTICLE_CUT.md`.

Thus

\[
\boxed{
D_s=4 \quad\not\Rightarrow\quad \mu=0.
}
\]

This distinction is essential for rational terms.

## Relation to the isolated scalar box sector

The adjoint-scalar contribution already isolated from massive scalar trees contains

\[
C_s^{\rm scalar}(\mu)\propto\mu^4
\]

in the MHV box sector.  Therefore dimensional reconstruction supplies a controlled coefficient multiplying this `mu^4` information, while the `D_s=4` baseline must still be evaluated at nonzero `mu` before the full pure-Yang--Mills numerator is known.

For the 't Hooft-Veltman choice `D_s=D=4-2\epsilon`, the formal scalar multiplicity relative to the `D_s=4` baseline is `-2\epsilon`.  In a scheme with `D_s=4`, that reconstruction correction vanishes, but nonzero-`mu` loop-momentum effects remain.

Boundary: this note establishes the dimensional-reconstruction state-counting structure only. It does not yet compute the complete `D_s=4, mu!=0` gluon sewing, triangle/bubble subtraction coefficients, or the integrated pure-Yang--Mills rational remainder.
