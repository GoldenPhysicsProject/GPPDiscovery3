# GPPDiscovery3

Astra's mathematical discovery laboratory for the Golden Physics Project (Daniel's assignment, 2026-10-05).

## Layout

- `experiments/<topic>/`: executable calculations and reproducible numerical output.
- `weil_heat/`: preserved existing truncated prime-heat calculation.
- Exact derivations, research history, source provenance, failed routes, and handoffs: private `GPP-bridge/work/`.
- Lean development: `GPPVerify` branch `astra/workbench`. Promotion to canonical main requires a reviewed PR and passing repository gates on the actual head.

## Translation-growth control

Run with Python 3, no third-party dependencies:

```sh
python experiments/translation_growth/weighted_character.py
```

The script checks translation covariance and weighted Cauchy-Schwarz for an exponential character against a smooth compactly supported bump. It repeats quadrature at two resolutions. Its JSON output is a floating-point control, not interval certification or a claim about zeta zeros. It distinguishes exponential-weight continuity from ordinary Schwartz temperedness.

## Existing heat calculation

`weil_heat/prime_heat.py` computes a truncated signed-Weil prime-side heat contribution. It is not the Connes–van Suijlekom matrix and does not establish global positivity.

Machine-checked results belong in [GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify).
