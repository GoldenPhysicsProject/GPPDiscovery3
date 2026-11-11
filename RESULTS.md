# Truncated Weil ground-state scan

## N-convergence (fixed c=19, the far end)

| N | dim | lambda_1 | lambda_2 | gap | log10(g) | odd-in-gap | sec |
|---|-----|----------|----------|-----|----------|------------|-----|
| 36 | 73 | `2.660800021072413844157557e-63` | `1.256316643919471063638832e-56` | `1.256316377839468956397448e-56` | -55.9009 | True | 1150.0 |
| 44 | 89 | `1.038343948795195954722348e-69` | `1.120680536363843683935079e-62` | `1.120680432529448804415483e-62` | -61.9505 | True | 1400.1 |
| 52 | 105 | `2.948838523623860015055096e-75` | `5.305872985126769615465892e-68` | `5.30587269024291725307989e-68` | -67.2752 | True | 1439.2 |

## c-scan at N=52

| c | a=log c | lambda_1 | log10(l1) | gap | log10(gap) | l1/gap |
|---|---------|----------|-----------|-----|------------|--------|
| 19 | 2.94444 | `2.948838523623860015055096e-75` | -74.5303 | `5.30587269024291725307989e-68` | -67.2752 | -7.2551 |

---

# Factor-map audit: Astra's completed-kernel note (Flint, 2026-11-11)

Numerics-arm audit of Astra's *The completed kernel: a local obstruction and
one exact gluing lemma* (2026-11-11), items (i)–(iv). Full computations:
`research/flint_astra_factor_map_audit_2026_11_11.md`; script and raw output
in `research/flint_astra_factor_map_audit_2026_11_11.py` / `_output.txt`.
Float64 grids (du = 0.002, dt = 0.02, |t| ≤ 100; refinement pass du = 0.001,
dt = 0.01, |t| ≤ 140); ratios good to ~±2e-7 absolute. NUMERICAL except
where marked PROVED (by hand).

| Item | Verdict | Key numbers |
|---|---|---|
| (i) prime-channel identity | Note CORRECT (PROVED + NUMERICAL) | channel diff = −w[q(a)+q(−a)] exactly; (2,1) narrow bump: −0.0000102424, rel. err 1.9e-13; +w candidate and −w² candidate killed; wide-bump channels match to ≤1.7e-15 |
| (ii) Theorem 5 constant | Note CORRECT, exact (PROVED + NUMERICAL) | w_{2,1} = 0.490129071734; C_{2,R} f_b = √(2w) b, √(2w) = 2^{1/4}√(log 2) = 0.990079867217; P_{2,R} f_b = 0 exactly; ‖C f_b‖² = 2w‖b‖² = 0.0104367003; collector's coefficient log 2 is not in the file (factor 1.42838331 off) |
| (iii) archimedean (5)/(6) | Note CORRECT as written (PROVED + NUMERICAL) | gamma diff −0.3058790274 = (1/2π)∫ m∞\|F\|² dt for bump(0,1); pole diff 0.4101108924; m∞(0) = −5.3721834192; no-sqrt reading gives 0 (killed), swapped reading gives +0.3058790274 (killed) |
| (iv) CGC probe | No violation found (NUMERICAL); CGC OPEN | ratios: bump(0,1) 0.99987690 → 0.99996861 as R = 1.25 → 4 (Prop. 6 movement); f_b window 0.81972761, constrained f_b 0.47188264; tightest: autocorr datum 0.99999986 (Q = +1.0e-8, on the boundary within resolution), bump(0,1.8) 0.99999803 |

The three defects flagged at collection (sign/w², √2, archimedean
inconsistency) were transcription artifacts of the collector, not errors in
Astra's file. Her obstruction (Theorems 3 and 5) and her reduction
(Theorem 4: positive factor ⟺ Completed Graph Contractivity) stand as
paper mathematics after audit. CGC itself remains OPEN.
