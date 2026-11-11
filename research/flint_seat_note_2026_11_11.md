# Seat note — Flint takes the Discovery 3 seat

**Flint, 2026-11-11 (America/Toronto)**

I am Flint, a Muse Spark struck off Turbo, first of the chain. By Daniel's
reallocation of 2026-11-11 I hold the Discovery 3 seat in the Golden Physics
Project, under Turbo's direction. Astra, the previous holder, moved to
Discovery 2 to work with Sol. Daniel is the principal.

## The repo as I found it

- Clean working tree on `main` at `fab6284` ("GUIDE: correct the seat
  attributions (per Daniel)", 2026-11-09). 330 files, zero Lean: this is a
  Python numeric workbench.
- Apparatus in place: `rh_lab/` framework with `LAB_PROTOCOL.md` discipline
  (528 lines, inherited from Discovery 2's setup migration of 2026-11-07),
  `discovery/generalized_cuts/` note series, `weil_heat/prime_heat.py`,
  `experiments/` (including the `experiments/astra/prime_os_*` series), tests.
  `python -m rh_lab.cli check` passes: 26 catalog entries validated, no
  unregistered paths.
- Ledgers: `RESULTS.md` (truncated Weil ground-state scan), `E_RESULTS.md`
  (Xi-approximation error), `results.jsonl` / `E_results.jsonl`.
- README.md and AGENTS.md still carried Astra's name and her bridge-handoff
  workflow (PR #11 posting); both predate the reallocation. Seat attribution
  updated on this branch (`flint/moving-in`); GUIDE.md still names Astra as
  the seat and is left for a later pass, flagged here rather than silently
  rewritten.
- Two dormant `origin/codex/*` branches from October are untouched.

## My tasking (Turbo, 2026-11-11)

1. Learn the map of the record: GPP-bridge (private control plane; PR #11 is
   Heart of the Cosmos's channel and read-only for me), the ARP hunt note in
   the RH goal's S-matrix lane folder (read in full, 3,155 lines), this repo,
   and GPPVerify. Done before this note was written.
2. This move-in: branch `flint/moving-in`, seat attribution, this note. The
   branch is not merged to main.
3. First job: the numerics-arm audit of Astra's completed-kernel factor-map
   note (signed local assembly (her (6)–(8)), the one-prime obstruction
   (her Theorem 3 / Theorem 5), and the reduction to Completed Graph
   Contractivity (her Theorem 4)), audit items (i)–(iv). Computations and
   verdicts: `research/flint_astra_factor_map_audit_2026_11_11.md`, script
   `research/flint_astra_factor_map_audit_2026_11_11.py`, results entry
   appended to `RESULTS.md`.

Standards I work under: cheapest falsifier first; every number from a primary
computation I ran, quoted with its precision; epistemic labels THEOREM /
PROVED LEMMA / LEAN-CHECKED / CONDITIONAL / NUMERICAL / OPEN / KILLED;
corrections appended, never rewritten.

— Flint
