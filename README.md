# GPPDiscovery3

Numeric and exploratory research workbench for the Golden Physics Project's shadow
holography framework — Astra’s discovery workbench alongside
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
[`discovery/arithmetic_holography/ARITHMETIC_HOLOGRAPHY_RIGOROUS.md`](https://github.com/GoldenPhysicsProject/GPPDiscovery/blob/main/discovery/arithmetic_holography/ARITHMETIC_HOLOGRAPHY_RIGOROUS.md).

**Nothing in this repo is proved.** It holds numeric evidence, hand derivations, and
literature checks. Once a result is solid enough to state as a theorem, it is formalized in
[GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify) (no `sorry`, no axiom
asserting an open claim, CI-verified). `AGENTS.md` documents that workflow and the
branch-hygiene rules.

## Running

Use Python 3.11 or later from the repository root:

```sh
python -m rh_lab.cli list
python -m rh_lab.cli check
python -m unittest discover -s tests -v
python -m rh_lab.cli smoke
```

The exact core uses the Python standard library. Optional experiments declare their dependencies in `rh_lab/catalog.json`; consult `LAB_PROTOCOL.md` before running numerical or interval checks.

## Layout and evidence

`experiments/`, `manuscripts/`, `research/`, and `weil_heat/` retain Astra's existing work. The runnable `rh_lab/` framework, tests, catalogued experiments, and CI setup are imported from Discovery2 with source provenance in `research/setup_migration_2026_11_07/PROVENANCE.md`. Imported result ledgers remain Discovery2 results, not new Astra measurements.

`RESULTS.md`, `E_RESULTS.md`, and the JSONL ledgers record results. `AGENTS.md` sets operating rules; `LAB_PROTOCOL.md` specifies assumptions, falsifiers, replay receipts, and promotion criteria. Detailed working notes are maintained in the private GPP bridge.

## Existing Astra calculations

```sh
python experiments/translation_growth/weighted_character.py
```

This floating-point translation-covariance control is not interval certification. `weil_heat/prime_heat.py` remains the truncated prime-heat calculation, not a proof of global positivity.

## Provenance

This front page follows Discovery1's README. Its full source and the previous Discovery3 README are preserved in `research/setup_migration_2026_11_07/`. Research-specific claims and scan plans from Discovery1 are retained there as source material, not attributed to Astra. No prior research files were deleted.
