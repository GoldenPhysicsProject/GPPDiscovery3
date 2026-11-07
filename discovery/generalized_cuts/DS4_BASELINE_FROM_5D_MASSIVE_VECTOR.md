# Nonzero-mu `D_s=4` baseline from the five-dimensional massive-vector embedding

Codex/GPT discovery track, 2026-08-25.

## Problem being closed

The dimensional-reconstruction identity already recorded is

\[
C^{(D_s)}(\mu)
=C^{(4)}(\mu)+(D_s-4)C^{(S)}(\mu),
\]

where `C^(S)` is one real adjoint-scalar state transported around the cut. The remaining opaque piece was the baseline `C^(4)(mu)` at nonzero transverse mass. It cannot be replaced by the strict massless four-dimensional helicity cut because `D_s=4` does not imply `mu=0`.

There is a cleaner integer-dimensional representation of that baseline.

## Five-dimensional embedding

At one loop, embed the loop momentum in five dimensions while keeping all external kinematics four-dimensional. A five-dimensional null loop momentum can be written

\[
L^A=(\ell^\mu,\mu),
\qquad L^2=0
\quad\Longleftrightarrow\quad
\ell^2=\mu^2.
\]

Thus its four-dimensional projection is a massive momentum of mass `mu`.

A massless gauge boson in `D_s=5` has

\[
D_s-2=3
\]

physical polarizations. Under the four-dimensional projection those three polarizations form precisely the polarization triplet of a four-dimensional massive vector. Denote the corresponding two-particle cut by

\[
C^{(V_m)}(\mu)\equiv C^{(D_s=5)}(\mu).
\]

The scalar dimensional-reconstruction law, now based at `D_s=5`, is therefore

\[
\boxed{
C^{(D_s)}(\mu)
=C^{(V_m)}(\mu)+(D_s-5)C^{(S)}(\mu).
}
\]

This is algebraically equivalent to the previous `D_s=4` form because setting `D_s=5` in that form gives

\[
C^{(V_m)}=C^{(4)}+C^{(S)}.
\]

Hence the formerly opaque baseline is

\[
\boxed{
C^{(4)}(\mu)=C^{(V_m)}(\mu)-C^{(S)}(\mu).
}
\]

The state count is the immediate consistency check:

\[
3-1=2=D_s-2\big|_{D_s=4}.
\]

For a general integer spin dimension the same count reads

\[
3+(D_s-5)=D_s-2.
\]

## Why this advances the MHV cut calculation

The scalar piece `C^(S)` is already explicit in the adjacent-MHV sector: it is built from two massive-scalar trees and contains the `mu^4` rational box numerator. Therefore the nonzero-`mu` `D_s=4` state sum can now be attacked by a concrete four-dimensional massive-vector sewing problem rather than by an abstract D-dimensional polarization tensor.

The immediate calculation becomes:

1. construct the two massive-vector tree amplitudes entering the `1^- 2^- | 3^+ 4^+` two-particle cut;
2. sum the three physical massive-vector polarizations;
3. subtract the already-known one-real-scalar cut;
4. only then combine with `(D_s-4) C^(S)` for the desired regularization scheme.

Equivalently,

\[
C^{(D_s)}=C^{(V_m)}+(D_s-5)C^{(S)}.
\]

This formulation is particularly natural for integer-dimensional generalized unitarity, where the loop-momentum dimension and spin-state dimension are kept distinct and amplitudes are reconstructed from integer-dimensional state sums.

## Scheme specializations

- `D_s=5`: `C^(D_s)=C^(V_m)`.
- `D_s=4`: `C^(4)=C^(V_m)-C^(S)`.
- Formal HV continuation `D_s=4-2 epsilon`:

\[
C^{HV}=C^{(V_m)}-(1+2\epsilon)C^{(S)}.
\]

This is equivalent to `C^(4)-2 epsilon C^(S)` after substituting `C^(4)=C^(V_m)-C^(S)`.

## Boundary

This note reduces the state-sum problem; it does not yet evaluate the massive-vector trees or their three-polarization cut. It also does not assert a generic box-only reduction. Triangle, bubble and subtraction sectors remain topology/helicity dependent. The important advance is that the previously unspecified `D_s=4, mu!=0` baseline now has a concrete four-dimensional massive-vector-minus-scalar realization.

References for the reconstruction framework: generalized D-dimensional unitarity treats the internal spin dimension `D_s` separately from loop-momentum dimension and uses integer-dimensional evaluations; at one loop a five-dimensional momentum embedding is sufficient. The scalar interpretation of extra-dimensional components and massive four-dimensional cut lines is standard in D-dimensional unitarity.
