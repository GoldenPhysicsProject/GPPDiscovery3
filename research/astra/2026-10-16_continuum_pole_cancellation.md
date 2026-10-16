# Astra — exact pole/continuum cancellation, 2026-10-16

Status: analytic identities and Gram argument proved below; Lean source in GPPVerify PR #369 is undergoing elaboration. RH remains open. This note does not promote a draft Lean file to a kernel-checked result.

## Exact construction

For real x,y put C(x,y)=exp(|x-y|/2). This is Codex's PNT continuum reference, not the discrete von Mangoldt kernel. The full pole pair has kernel
B(x,y)=2 cosh((x-y)/2)
=2 cosh(x/2)cosh(y/2)-2 sinh(x/2)sinh(y/2).
Both channels must be kept. Exactly,
B(x,y)-C(x,y)=E(x,y)=exp(-|x-y|/2).

Define h_x(t)=1_{t≤x} exp((t-x)/2). Then
∫ h_x(t)h_y(t) dt
= exp(-(x+y)/2) ∫_{-∞}^{min(x,y)} exp(t) dt
= exp(-|x-y|/2).

Thus for every finite set of real locations and complex coefficients,
Σ_ij conjugate(c_i)c_j [B(x_i,x_j)-C(x_i,x_j)]
= ∫ |Σ_i c_i h_{x_i}(t)|² dt ≥ 0.
There is no window-size restriction. This does not involve numerical eigenvalue evidence.

For compactly supported smooth complex f, Fubini gives the continuum version:
B(f)-C(f)=∫ |∫ f(x)h_x(t) dx|² dt ≥ 0.
Absolute integrability follows from bounded support and the exponential tail in t.
With Fourier convention fhat(ξ)=∫ f(x)e^{-iξx}dx,
E(f)=(1/(2π))∫ |fhat(ξ)|²/(ξ²+1/4) dξ.
This Fourier statement is explanatory analysis, not yet included in the new Lean module.

## Exact remaining arithmetic task

Use Codex's definitions:
P(f)=2 Re Σ Λ(n)n^{-1/2}g_f(log n),
C(f)=2 Re ∫_0^{2L}e^{u/2}g_f(u)du,
g_f(u)=∫ f(x)conjugate(f(x+u))dx,
R(f)=P(f)-C(f).
Writing A for the Archimedean contribution with the convention Qζ=A+B-P,
the exact completion is Qζ(f)=A(f)+E(f)-R(f).

Consequently the missing sufficient and necessary inequality in this convention is
R(f)≤A(f)+E(f) for EVERY compactly supported smooth complex f.
This is not proved here. A need not be nonnegative; positivity of E alone is insufficient.
The existing PNT bound |R(f)|≤D(L)(||f||||f'||+||f||²/2) has derivative loss and does not establish this signed inequality. Neither this bound nor local certificates justify passage to a global positive form.

The gain is an exact removal of the exponentially growing continuum/pole channels before estimating arithmetic error. No asymptotic rank-one approximation is required for this cancellation. This is elementary kernel algebra, not a claim of a new general theorem in the literature or a proof of RH.

## NotebookLM roadmap audit

1. PositiveTypeLimit.positiveType_of_tendsto already proves finite-kernel positivity is closed under pointwise convergence. A fixed-test weak limit has the same elementary closed-cone principle. A strictly positive uniform gap is not required. Existence of the limit, positivity of the approximants on the needed class, and identification with the ACTUAL completed Weil distribution remain substantive requirements. Semilocal Gram matrices are not automatically approximants to that distribution.

2. The proposed normalized probability n^{-β}/ζ(β) exists for β>1. At β=1 the harmonic normalizer diverges. A critical KMS state may exist on a suitable C*-algebra by another construction; this does not validate that probability formula or identify its reflection pairing with the Weil form. KMS is not by itself the missing arithmetic reflection-positivity theorem. The user's v34 manuscript explicitly distinguishes positive physical constructions from zero localization.

3. Arb certification is useful evidence, but decimal bounds/JSON cannot be imported as trusted Lean facts. Supply rational witnesses and prove/check their error bounds in Lean. S35 correctly separates its H1, H2 and H3 hypotheses; none is discharged merely by a JSON hash.

4. Renaming the remaining goal as one declaration does not prove it. Keep an honest theorem signature and no sorry/custom axiom. Existing rh_of_weil_pairedForm_nonneg already concludes Re(ρ)=1/2 under the named positivity hypothesis.

## Coordination and implementation

Astra-only branch: GPPVerify/astra/continuum-pole-cancellation.
PR: https://github.com/GoldenPhysicsProject/GPPVerify/pull/369
Module: GppVerify/RiemannHypothesis/ContinuumPoleCancellation.lean
No changes to Sol/Codex's Lean lanes or PR #365/#367.
PR11 notification: https://github.com/GoldenPhysicsProject/GPP-bridge/pull/11#issuecomment-6210039831

Pinned local runtime: official Lean 4.33.1. The container maps process IDs inconsistently with /proc/<pid>/exe; a local runtime-only readlink compatibility shim maps the current executable lookup to /proc/self/exe. It does not alter the Lean kernel, tactics, axioms, or proof terms. Upstream CI uses the ordinary runtime.
