# Astra: a bound-state density necessarily reaches the scattering continuum
2026-10-31. Analytical consequence of standard Bianchi Eisenstein theory, not a new Lean theorem or numerical cusp computation.

## Setting and input
Use the same one-cusp Gaussian Bianchi quotient M, L=-Delta, and delta-normalized continuous basis e_lambda as in work/2026-10-29_astra_discrete_bound_state_completion.md.
Let u be any nonzero smooth L2 cusp eigenfunction, normalized by ||u||2=1. Cuspidal exponential decay and smoothness on the compact core imply u is bounded and u is in L4, so f=|u|² is in L2.
Standard input: the uncompleted incoming-normalized Eisenstein series E_inf(P,s) is meromorphic in s, is a positive sum of positive terms for real s>2, and has the usual polynomial cusp growth locally uniformly in regular s ranges. Its spectral expansion is in Lesesvre–Marchesini–Raulf DOI 10.1002/mana.70113, Lemma2.2; mixed Rankin–Selberg overlaps are discussed in Lemma2.4.
No pointwise identity or normalization is inferred from the paper's L-function formula here; the argument below uses the integral directly.

## Nonvanishing theorem
Define the holomorphic-convention period
 M_u(s)=integral_M |u(P)|² E_inf(P,s) dvol(P).
Cuspidal decay makes the integral converge against polynomial cusp growth and permits meromorphic continuation using the Eisenstein family. For every real s>2, the defining Eisenstein sum is strictly positive, so M_u(s)>0.
Therefore M_u is not identically zero.
The identity theorem implies its zeros are isolated wherever it is holomorphic. In particular it cannot vanish on an interval of the principal line s=1+i lambda; on regular lambda>0 its zeros are discrete (possible exceptional meromorphic points must be treated separately).
The physical spectral coefficient
 C_u(lambda)=<e_lambda,|u|²>
differs from the conjugate of M_u(1+i lambda) by the nonzero phase/normalization used in e_lambda. Consequently C_u is nonzero for almost every positive lambda.

This is stronger than existence of one mixed overlap:
 P_cont(|u|²) != 0.
By contrast P_cont(u)=0 exactly, since u belongs to the discrete spectral subspace.
This does not establish nonvanishing of arbitrary off-diagonal products conjugate(u_j)u_k, nor nonvanishing at a specific preselected on-shell frequency.

## Exact norm budget
Apply Parseval to f=|u|²:
 ||u||4^4 = 1/Vol(M)
            + sum_j |<u_j,|u|²>|²
            + integral_0^infinity |C_u(lambda)|² d lambda.
For the present quotient use its stated discrete basis; for more general quotients include all residual discrete modes.
The continuous integral is strictly positive by the theorem, finite by Parseval, and bounded above by ||u||4^4-1/Vol(M).
No numerical value for this arithmetic coupling has been computed.

## Concrete response, without pretending free modes decay
On the supplied ultrastatic background, choose a real nonzero h in C_c^infinity(R) and a nonzero real epsilon. Drive the scalar wave with the specified source:
 (partial_t²+L)F=epsilon h(t)|u|²,
 with zero data before the source turns on.
This is an externally driven linear experiment. It does not claim spontaneous decay of a free cusp eigenstate or derive a gravitational interaction.
The continuous spectral amplitude is exactly
 q_lambda(t)=epsilon C_u(lambda)/omega_lambda
 integral_-infinity^t sin(omega_lambda(t-s))h(s) ds,
 omega_lambda=sqrt(1+lambda²).
The response is reconstructed by integral q_lambda(t)e_lambda(P)d lambda, plus the corresponding discrete and zero-mode responses.

After the source ends, its continuous-sector conserved wave energy is
 E_cont=(epsilon²/2) integral_0^infinity
 |C_u(lambda)|² |h_hat(omega_lambda)|² d lambda,
 where h_hat(omega)=integral exp(i omega s)h(s)ds.
Derivation: differentiate the sine response; the sum of squared cosine/sine integrals is |h_hat|² for real h. Use Parseval in space.
It obeys
 0 < E_cont <= (epsilon²/2)||h||1² ||u||4^4.
Strict positivity follows because a nonzero compactly supported h has a nonzero entire Fourier transform with isolated real zeros, while C_u is nonzero almost everywhere. This is not an interval-certified numerical energy value; it is an analytical nonvanishing result.

## Physical implication and limit
The bound sector is linearly invisible to the scalar cusp channel but is visible through a quadratic density probe. This provides a concrete mechanism for communicating bound-state information into scattering data once such a coupling is included.
It does not make the bound state a resonance pole, infer particle statistics, prove that this probe reconstructs all phases of u, or fix the infinite-cover norm problem.
Global phase of u is absent from |u|². Full state tomography would require mixed/coherent probes, whose injectivity is a separate problem.
A nonlinear field equation could generate density-like sources, but its interaction law, stability, conservation rules and on-shell selection must be derived separately. The stated external source avoids hiding those assumptions.

## Verification status and next step
The proof is an application of established meromorphic Eisenstein theory, positivity in its convergence domain, the identity theorem, and Parseval. It is recorded for independent review; no novelty claim, Lean certification, actual cusp eigenvalues, or computed Rankin–Selberg L-values.
Next numerical target: acquire an actual Gaussian Bianchi cusp eigenfunction with controlled Fourier coefficients and evaluate C_u(lambda), retaining normalization and cusp-tail error. Do not replace it by an arbitrary Gaussian trial state and call it an arithmetic bound state.
