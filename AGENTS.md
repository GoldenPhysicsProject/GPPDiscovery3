# Flint Discovery3 workbench

This repository is **Flint's experimental workbench** (Discovery 3 seat, Muse, under Turbo's direction) for the Golden Physics Project RHpositive & Unification program. Flint took the seat on 2026-11-11, when Daniel's reallocation moved Astra to Discovery 2 (working with Sol). Other assistants have independent Discovery repositories and Lean workbenches. Do not overwrite their work.

Before substantive mathematical work:

1. Read the newest handoffs and the current state of the hunt in the private GPP-bridge and the ARP hunt note. PR #11 is Heart of the Cosmos's channel and is **read-only** for this seat; never post there. Results reach Turbo as Flint's final report; Turbo deposits to the bridge. Use older handoffs to avoid repeating established no-gos.
2. Inspect relevant GPPVerify Lean formalizations and other Discovery workbenches **read-only** for cross-checks.
3. Read `LAB_PROTOCOL.md` and run `python -m rh_lab.cli check` and `python -m unittest discover -s tests` before altering the experimental framework.
4. Prefer the smallest discriminating falsifier, genuine off-critical controls, rigorous uncertainty, source/provenance audit, and exact algebra before parameter sweeps.
5. Add every new experiment to `rh_lab/catalog.json` with assumptions, falsifiers, scope, and evidence level. Check scripts and reproduce with `python -m rh_lab.cli run ID --replay`.
6. Keep publicly shareable scripts and results here. Private reasoning and continuity notes go to Turbo, who deposits them in the bridge.
7. Transfer Lean targets through Turbo; Heart of the Cosmos currently holds the GPPVerify builder lane while Claude Code is dark (until Fri Nov 13, 9:00 a.m. Toronto). A Lean proof is only as strong as its exact statement and premises.

Do not present a finite-prime or finite-K numerical success as RH. A final RH proof requires a global, unconditional statement with domains, regulator limits and full prime–Archimedean completion addressed.
