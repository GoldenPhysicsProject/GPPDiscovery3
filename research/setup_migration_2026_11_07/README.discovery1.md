# GPPDiscovery

Numeric and exploratory research workbench for the Golden Physics Project's shadow
holography framework — the standing discovery repo behind
[GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify), which carries the
formalized, Lean 4-checked results.

## The framework, briefly

**Shadow holography.** One involution appears in two places. It is the celestial shadow map
`Δ ↦ 2 − Δ` on conformal dimensions, and the zeta functional equation `s ↦ 1 − s`. Under
`Δ = 2s` they coincide, and the reflection `τ(s) = 1 − s̄` fixes exactly the critical line
(celestial `Re Δ = 1`). The framework, developed in Daniel Toupin's *On the Nature of
Nature*, reads this as a nested holographic dictionary. Each level is a Mellin or
Pontryagin duality, so each is exact:

- **Line ↔ sphere.** The multiplicative line `ℝ₊`, split at `x = 1`, is the boundary. The
  spectral sphere, split by the critical line (the equator), is its hologram. The inversion
  `x ↦ 1/x` becomes `s ↦ 1 − s̄`.
- **Integers ↔ functional equation.** The integers are the Fourier modes of a circle, and
  Poisson self-duality becomes the functional equation.
- **Primes ↔ torus.** The primes span a compact torus `∏ S¹`, and the critical line embeds
  in it densely: a hologram inside a hologram.

In this language the functional equation says the hologram is *symmetric*. RH is one open
question inside the framework: whether the hologram is also *positive* (reflection-positive
across the equator). Celestial scattering meets the same equator in its scale sector. That
is a structural parallel, not a claim that `ξ` is a physical correlator. The full write-up,
with what is a theorem, what is open, and what has been killed, is
[`discovery/arithmetic_holography/ARITHMETIC_HOLOGRAPHY_RIGOROUS.md`](discovery/arithmetic_holography/ARITHMETIC_HOLOGRAPHY_RIGOROUS.md).

**Nothing in this repo is proved.** It holds numeric evidence, hand derivations, and
literature checks. Once a result is solid enough to state as a theorem, it is formalized in
[GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify) (no `sorry`, no axiom
asserting an open claim, CI-verified). `CLAUDE.md` documents that workflow and the
branch-hygiene rules.

## Active threads

- **`discovery/arithmetic_holography/`**: making the shadow-holography dictionary rigorous.
  This covers the three levels above, the Osterwalder–Schrader form of the one open
  positivity statement, and the constraints any route must respect (Beurling /
  Davenport–Heilbronn, de Bruijn–Newman, the Fisher no-go). It includes a numerical
  falsifier (`riemann_split_hb.py`) that killed one proposed hemisphere-dominance route.
- **`weil_decay/`** (root scripts + `discovery/weil_decay/`) — the truncated Weil
  quadratic form. Connes–van Suijlekom and Connes–Consani–Moscovici build, for a prime
  cutoff `c` and band `N`, a finite Galerkin matrix `Q(c)` whose zeros provably sit on
  the critical line for every finite `c`; convergence as `N → ∞` is the open question
  this thread measures. See "The Weil-decay question" below for the current numbers.
- **`discovery/celestial_box/`** — extending the audited scalar-box cut/dispersion
  construction toward pure Einstein gravity; tracking exactly where D-dimensional
  unitarity subtleties (rational terms invisible to 4D cuts) start to matter.
- **`discovery/wiener_hopf/`** — the exact bridge between the principal-series cut
  weight, the Wiener–Hopf Fourier window, and the conical prefactor at
  `s = 1/2 + it`.
- **`discovery/positive_reals_cft/`** — isolating the precise representation-theoretic
  structure behind treating `(R⁺, ×)` as a one-dimensional principal-series system.
- **`discovery/number_thermodynamics/`** — the canonical Gibbs distribution
  `P_β(n) = n⁻β/ζ(β)` on the positive integers and its thermodynamic reading.

## Thread detail: the Weil-decay scan

`λ_min(c)`, the smallest even-sector eigenvalue of `Q(c)`, is non-negative for every
finite `c`. (Weil positivity, and hence RH, is this holding in the limit. This is one numerical
probe of the positivity question, not a route the project is built around.) We are not testing
whether it's positive — we're measuring **how fast it decays**:

> Is `log λ_min` linear in `log c`, and if so, what is the constant?

Measured so far (unconverged, N=20, T=100): `3.67e-27` at c=7, `8.09e-18` at c=5. A
two-point natural-log slope of about `-64`, against `-2γ₁ = -28.27` and
`-4γ₁ = -56.55`. If the constant is arithmetic, this is a quantity that *sees the
zeros*, which is more than can be said for most reformulations of RH.

**N-convergence comes first.** Going from dim 33 to dim 41 at c=7 moved `λ_min` from
`2.4e-25` to `3.67e-27`, a factor of 65 — N=20 is not converged and any slope fitted to
it is meaningless. The scan does c=7 at N = 20, 28, 36 before anything else.

**Caveat on small c.** At c=3 (N=20, T=100) we saw `λ_min = -0.332` — almost certainly
marginal basis resolution, per the reference implementation's own account at c=23, 29.
The c=29 case has since been resolved at higher precision: stably positive
(`≈1.59e-62`, 120-digit residual `1.64e-121`) at 90–150 digits, though N and T are still
unconverged there too, so this doesn't yet rescue a c-slope fit. See
[`discovery/weil_decay/C29_PRECISION_AUDIT.md`](discovery/weil_decay/C29_PRECISION_AUDIT.md).

## Running

`python point.py --c 7 --N 36 --T 300 --dps 90` for one point, or push / use
`workflow_dispatch` to run the whole matrix in parallel on Actions. Results land in
`results.jsonl` and `RESULTS.md`, committed back automatically.

## Status

Discovery only, across every thread above. See `CLAUDE.md` for the process rules
(branch hygiene, when and how a result graduates to GPPVerify).
