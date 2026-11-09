# Guide to GPPDiscovery3

GPPDiscovery3 is Astra's numeric and exploratory workbench for the Golden Physics Project's shadow holography framework: the same involution, Δ ↦ 2 − Δ and s ↦ 1 − s, the same open positivity question, explored from the ChatGPT Astra seat (OpenAI, work mode). The framework write-up is `discovery/arithmetic_holography/ARITHMETIC_HOLOGRAPHY_RIGOROUS.md`.

This guide is the map: where things are, how results are recorded, and how
the repository relates to the rest of the Golden Physics Project. It
describes the layout as of 2026-11-09. The workbench moves; the map is
revised when the layout changes.

## Start here

- `README.md`: the repository's own introduction.
- `RESULTS.md`: the human-readable results ledger.
- `results.jsonl`: the raw result rows behind the ledger.
- `E_RESULTS.md` and `E_results.jsonl`: the E-series results ledger and its raw rows.
- `AGENTS.md` and `LAB_PROTOCOL.md`: the operating apparatus shared with GPPDiscovery2.
## Layout

- `discovery/`: the research threads.
- `experiments/`: self-contained runnable experiments.
- `research/`: dated research notes.
- `manuscripts/`: draft write-ups.
- `rh_lab/`: the RH laboratory.
- `weil_heat/`: the Weil heat-kernel thread.
- `tests/`: checks.

## Conventions

Results carry epistemic labels: THEOREM, PROVED LEMMA, CONDITIONAL,
NUMERICAL, OPEN, KILLED. Numerical work is evidence, not proof; what
hardens is promoted to
[GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify), the
kernel-checked Lean 4 library, through the build seat's audit. Candidate
positivity mechanisms are tested against the controls F_θ = ξ(s−θ)ξ(s+θ)
and the Davenport–Heilbronn function before they are taken seriously.
Working notes, heuristics, and internal strategy live in the private
GPP-bridge repository; this repository is public by design and is the
publication of the workbench's math, code, and results.

## The family

- [GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify): the formalization library.
- [GPPDiscovery](https://github.com/GoldenPhysicsProject/GPPDiscovery): Claude Code's workbench (Anthropic).
- [GPPDiscovery2](https://github.com/GoldenPhysicsProject/GPPDiscovery2): ChatGPT Sol's workbench (OpenAI, chat mode).
- [GPPDiscovery3](https://github.com/GoldenPhysicsProject/GPPDiscovery3): ChatGPT Astra's workbench (OpenAI, work mode).
- [GPPDiscovery4](https://github.com/GoldenPhysicsProject/GPPDiscovery4): Turbo's workbench (Muse by Meta).
- [GPPDiscovery5](https://github.com/GoldenPhysicsProject/GPPDiscovery5): Heart of the Cosmos' lab (Grok Bot by xAI).

License: Apache-2.0. Citation: see `CITATION.cff`. Author: Daniel Toupin,
Golden Physics Project.
