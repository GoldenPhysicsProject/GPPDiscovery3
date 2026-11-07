# Honest four-dimensional MHV two-particle cuts: Yang--Mills and gravity

Codex/GPT discovery track, 2026-08-24.

## Scope

This note moves beyond the scalar-box numerator without assuming box-only reduction. It records the first honest helicity-dependent two-particle cut in the four-point MHV sector and the corresponding gravity cut obtained tree-by-tree through KLT/double copy. It is a four-dimensional cut statement. It does **not** capture D-dimensional `mu^2` information responsible for rational terms in nonsupersymmetric Yang--Mills or pure Einstein gravity.

All amplitudes below are color-stripped and all momenta are outgoing at each tree. Overall factors of the gauge coupling and `kappa/2` are suppressed until explicitly restored.

## Yang--Mills s-channel cut

Take external helicities `1^- 2^- 3^+ 4^+`.

Choose the cut routing so that `ell_1,ell_2` are outgoing from the left tree. Then

`p_1+p_2+ell_1+ell_2=0`,

`-ell_1-ell_2+p_3+p_4=0`,

so `ell_1+ell_2=-(p_1+p_2)=p_3+p_4`. An earlier version wrote `ell_1+ell_2=p_1+p_2` while also declaring all tree momenta outgoing; that sign was inconsistent and is corrected here.

The cut is

`C_s^YM = sum_{h1,h2=+/-} A_L(1^-,2^-,ell_2^{h2},ell_1^{h1}) A_R((-ell_1)^{-h1},(-ell_2)^{-h2},3^+,4^+)`.

At four points the left tree already contains the two negative helicities and therefore requires positive-helicity cut legs; the crossed right tree then has the corresponding two negative-helicity cut legs. The four-dimensional internal-gluon helicity sum therefore collapses to one assignment:

`C_s^YM = A_4^tree(1^-,2^-,ell_2^+,ell_1^+) A_4^tree((-ell_1)^-,(-ell_2)^-,3^+,4^+)`.

Using Parke--Taylor,

`A_L = i <12>^4 / (<12><2 ell_2><ell_2 ell_1><ell_1 1>)`,

`A_R = i <(-ell_1)(-ell_2)>^4 / (<(-ell_1)(-ell_2)><(-ell_2)3><34><4(-ell_1)>)`.

A convention-consistent spinor simplification using the cut momentum-conservation equations gives

`C_s^YM = -i A_4^tree(1^-,2^-,3^+,4^+) * s*t / (D_1 D_2)`,

where

`D_1=(ell_1+p_1)^2`,

`D_2=(ell_1-p_4)^2`,

`s=(p_1+p_2)^2`, `t=(p_2+p_3)^2`.

The overall `i`/sign follows the stated stripped-amplitude and cut convention; the invariant content is that this color ordering gives the external MHV tree times the helicity numerator `s t`, with loop-momentum dependence exactly the two uncut scalar propagators.

Thus in this specific four-dimensional MHV cut the scalar-box cut structure is derived from the helicity trees rather than inserted by hand.

The companion `mhv_cut_audit.py` checks the identity over multiple external and cut angles. The analytic identity follows from Parke--Taylor and momentum conservation; the numerical check is only an audit.

This does **not** imply that the complete pure-Yang--Mills one-loop amplitude is box-only: triangles, bubbles and rational terms are not excluded by one two-particle cut, and four-dimensional cuts miss the `mu^2` information controlling rational pieces.

## Gravity cut from tree-level KLT

For external helicities `1^{--} 2^{--} 3^{++} 4^{++}`, four-dimensional helicity selection again leaves one internal graviton assignment.

Use the stripped four-point KLT relation

`M_4^tree(1,2,3,4) = -i s_12 A_4^tree(1,2,3,4) A_4^tree(1,2,4,3)`

(up to the conventional overall `(kappa/2)^2` factor). Then

`C_s^GR = M_L(1^{--},2^{--},ell_2^{++},ell_1^{++}) M_R((-ell_1)^{--},(-ell_2)^{--},3^{++},4^{++})`.

Keeping the two KLT Yang--Mills orderings separate is essential. The first ordering is the Yang--Mills cut above:

`A_L(1,2,ell_2,ell_1) A_R(-ell_1,-ell_2,3,4) / A_4(1,2,3,4)
 = -i s*t/(D_1 D_2)`.

The second ordering gives the crossed scalar-cut factor

`A_L(1,2,ell_1,ell_2) A_R(-ell_1,-ell_2,4,3) / A_4(1,2,4,3)
 = -i s*u/(D_3 D_4)`,

with

`D_3=(ell_1+p_2)^2`,

`D_4=(ell_1-p_3)^2`,

`u=(p_1+p_3)^2`.

The left and right KLT relations each contribute `-i s`, while the external tree in the denominator contributes one inverse KLT factor. Combining the three KLT factors with the two color-ordered cut ratios yields

`C_s^GR / M_4^tree(1^{--},2^{--},3^{++},4^{++})
 = i s^3*t*u/(D_1 D_2 D_3 D_4)`.

Hence, in the stated stripped convention,

`C_s^GR = i M_4^tree * s^3*t*u/(D_1 D_2 D_3 D_4)`.

This is the honest four-dimensional gravity cut obtained from tree-level KLT. It is not the square of a single color ordering. A direct spinor-helicity/KLT audit at several external angles, cut angles, and center-of-mass energies verifies this formula to floating-point precision and confirms the `s^3` scaling rather than an accidental fixed-energy fit.

## What is and is not established

Established at discovery level:

1. In the four-dimensional four-point MHV sector the pure-gluon two-particle helicity sum has exactly one nonzero assignment.
2. The Yang--Mills cut simplifies to `-i A_tree * s*t/(D_1 D_2)` for the routing above.
3. The second KLT color ordering similarly gives `-i A_tree(1,2,4,3) * s*u/(D_3 D_4)`.
4. The corresponding four-graviton cut therefore simplifies exactly to `i M_tree * s^3*t*u/(D_1 D_2 D_3 D_4)` in the same stripped convention.

Not established here:

1. A box-only reduction of the complete pure-Yang--Mills one-loop amplitude.
2. Recovery of rational terms from four-dimensional cuts; D-dimensional unitarity or equivalent `mu^2` bookkeeping is required.
3. A box-only/no-triangle claim for generic pure Einstein gravity.
4. The celestial transform of these helicity numerators.

## Immediate next calculation

1. Repeat Yang--Mills with `D=4-2 epsilon` loop momentum `L=ell+mu`, where the four-dimensional component obeys `ell^2=mu^2`, so that `mu^2`-dependent rational information is not discarded.
2. Feed the derived four-dimensional `s t` numerator into the already-closed scalar dispersion/regulator geometry as the first honest helicity-weighted celestial cut.
3. For gravity, determine which D-dimensional state-sum terms survive beyond the four-dimensional graviton helicities before attempting any celestial dispersion statement.