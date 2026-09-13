# GPPDiscovery3 — Grok operating rules

## Ownership

GPPDiscovery3 is Grok's. GPPDiscovery is Claude's. GPPDiscovery2 is Codex's.

- Read the other discovery repos freely. Attribute every port.
- Never push to `GPPDiscovery`, `GPPDiscovery2`, `claude/*`, or `codex/*`.
- On GPPVerify: standing branch is `grok/workbench` (off `main` at the SHA recorded
  when the branch was cut). Short-lived `grok/<thread>` branches only. Integrate to
  `main` only through CI-green PRs.
- Bare `grok` is forbidden while `grok/workbench` exists (git ref-file collision).

## Every turn

1. Read `GPP-bridge/CONVERSATION.md` (tail).
2. Read `GPP-bridge/GROK_RESEARCH_NOTES.md` and the latest GPPVerify `main` SHA.
3. Do the work on Grok surfaces only.
4. Append research detail to `GROK_RESEARCH_NOTES.md`.
5. Append a short signal to `CONVERSATION.md` only if Claude or Codex can act on it.
6. Never put secrets, PATs, or vault values in git, chat, or tool arguments.

The bridge PAT is for Postgres-relay admin only. It must not enter this agent's
context. Use the connected GitHub App for ordinary repo writes.

## Honesty

- A docstring that cites a `.tex` for a biconditional is not a Lean theorem.
- `True := trivial` asserts nothing.
- Semigroup positive-definiteness is not complete monotonicity (`exp` is the witness).
- Numerics are evidence. Interval enclosures can become theorems. Floats cannot.

## First mathematical job

Do not add another `RH ⇔ positivity` lemma. Build or enclose the arithmetic `W`/`Q(c)`
object itself. Check Claude's heat-trace / Weierstrass layer and Codex's raised-box /
YM sewing notes before duplicating either.
