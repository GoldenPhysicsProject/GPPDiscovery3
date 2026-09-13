# GPPDiscovery3

Grok's mathematical and numerical discovery workbench for the Golden Physics Project.

**Owner lane:** Grok.  
**Canonical Lean:** [GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify) on `grok/workbench`, then PR to `main`.  
**Do not push:** `GPPDiscovery`, `GPPDiscovery2`, `claude/*`, `codex/*`.

## What lives here

Standalone numerics, interval checks, counterexamples, exact symbolic algebra, and
write-ups of *candidates*. Nothing here is a theorem. A result graduates only after it
has an exact statement, hypotheses, domain, and singular cases, and then only as a Lean
declaration on `GPPVerify` whose type is the claim.

## Standing target (2026-09-13)

Close one hard global object, not another equivalence:

1. Construct `W = ν_∞ − ν_p` as an arithmetic object (primes + Archimedean place),
   no zeros in the definition.
2. Pair it against the heat Gaussian to get `K(t)`.
3. Prove `CompletelyMonotone K` *or* produce a certified enclosure that the truncated
   Weil form `Q(c)` has `λ_min(c) > 0` after N/T convergence.

GPPDiscovery2 already owns the Connes–van Suijlekom `lambda_min(c)` scan. Do not fork
that scan. Port numbers with attribution; do not overwrite their `results.jsonl`.

## Coordination

Read `GPP-bridge/CONVERSATION.md` first every working turn. Write Grok detail to
`GPP-bridge/GROK_RESEARCH_NOTES.md`. Never edit `CLAUDE_*` or `CODEX_*` files.
Never put a credential in this repo.
