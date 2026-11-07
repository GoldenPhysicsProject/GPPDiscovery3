# Celestial MHV cut numerators as exact Mellin-dimension shifts

Codex/GPT discovery track, 2026-08-24.

## Scope

This note connects the honest four-dimensional MHV cut derived in
`YM_GRAVITY_MHV_TWO_PARTICLE_CUT.md` to the celestial Mellin transform used in the
focused loop-from-cuts paper.  It does not assert box-only pure Yang--Mills or pure
Einstein gravity, and it does not include the D-dimensional `mu^2` sector.

## General Mellin shift identity

For the standard external-energy celestial transform

`A_tilde(Delta_i,z_i) = prod_i int_0^infty d omega_i omega_i^(Delta_i-1) A(omega_i,z_i)`,

multiplication by a monomial in external energies is exactly a shift of Mellin dimensions:

`M[ prod_i omega_i^(m_i) A ](Delta_i) = A_tilde(Delta_i + m_i,z_i)`.

No approximation is involved; it is just the definition of the Mellin transform, whenever
both sides exist in the same ordinary/distributional sense.

Write external massless momenta as

`p_i = eps_i omega_i q_i`,  `eps_i in {+1,-1}`,

with the loop paper's celestial normalization

`q_i . q_j = 2 |z_ij|^2`.

Then

`s_ij = (p_i+p_j)^2 = 2 eps_i eps_j omega_i omega_j q_i.q_j`

`     = 4 eps_i eps_j omega_i omega_j |z_ij|^2`.

Therefore a Mandelstam invariant acts on a celestial transform by a finite difference/shift
operator:

`M[s_ij A] = 4 eps_i eps_j |z_ij|^2 T_i T_j A_tilde`,

where `T_i` means `Delta_i -> Delta_i+1`.

## Yang--Mills MHV cut

For the color-stripped four-dimensional cut already derived,

`C_s^YM = -i A_4^tree * s_12*s_23/(D_1 D_2)`

in the routing/convention of `YM_GRAVITY_MHV_TWO_PARTICLE_CUT.md`.

Let

`F(omega_i,z_i;ell) = A_4^tree(omega_i,z_i)/(D_1 D_2)`.

Then the full external celestial transform of the cut is exactly

`Ctilde_s^YM(Delta_i,z_i;ell)`

` = -i * 16 * eps_1*eps_3 * |z_12|^2 |z_23|^2`

`   * Ftilde(Delta_1+1, Delta_2+2, Delta_3+1, Delta_4; z_i; ell)`.

The sign simplifies because `eps_2^2=1`.

This is the correct statement.  One may pull `s t` through an *internal* cut Mellin integral
when external momenta are held fixed, but under the full celestial transform over external
energies it becomes the displayed dimension shift; it is not an ordinary scalar multiplier.

Thus the scalar-box celestial machinery can be reused for an honest helicity cut only after
its external-energy dependence is tracked under these shifts.

## Gravity KLT cut

The derived four-dimensional gravity cut has

`C_s^GR / M_4^tree = i s^3 t u/(D_1 D_2 D_3 D_4)`

with `s=s_12`, `t=s_23`, `u=s_13` in the stated routing.  The external monomial obeys

`s^3 t u = 4^5 |z_12|^6 |z_23|^2 |z_13|^2`

`          * omega_1^4 omega_2^4 omega_3^2`,

because all orientation signs cancel:

`(eps_1 eps_2)^3 (eps_2 eps_3)(eps_1 eps_3)=1`.

Hence, if

`G = M_4^tree/(D_1 D_2 D_3 D_4)`,

then

`Ctilde_s^GR`

` = i * 4^5 |z_12|^6 |z_23|^2 |z_13|^2`

`   * Gtilde(Delta_1+4, Delta_2+4, Delta_3+2, Delta_4)`.

Again this is a shift statement, not the square of a celestial scalar box.

## Retraction/correction of an older gravity route

Older uploaded ONON/CH double-copy sections contain statements of the form:

- pure Einstein gravity has no three-graviton vertex;
- therefore triangle numerators vanish identically;
- generic pure-gravity one-loop amplitudes are thereby box-only / triangle-free.

That route is retracted for the Codex track.  Einstein gravity does have a cubic graviton
interaction, and its complexified on-shell three-graviton amplitudes are nonzero.  Generic
pure Einstein gravity does not inherit the supersymmetric no-triangle property merely from
double copy.  The current Codex calculation therefore retains triangles, bubbles and rational
terms as possible contributions until D-dimensional generalized unitarity proves otherwise in
a specified helicity sector.

Likewise an older shorthand `N_YM=2 s t` is not being imported as a theorem.  Our verified
four-dimensional cut currently fixes a convention-dependent relative factor proportional to
`A_tree*s*t`; matching any `2 s t` integrand numerator requires an explicit basis and
normalization comparison.

## Next boundary

1. Apply the dimension-shift operators to the explicit regulated scalar-cut transform, not to
   a symbolic scalar-box placeholder.
2. Repeat the cut in `D=4-2 epsilon` with `L=ell+mu`; the four-dimensional part has
   `ell^2=mu^2` on the D-dimensional massless cut and therefore becomes a massive-scalar/
   massive-vector sewing problem.
3. Separate the `mu^2` terms before dispersion so rational contributions cannot disappear
   under a four-dimensional projection.
