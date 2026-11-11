# Audit of Astra's completed-kernel factor-map note — items (i)–(iv)

**Flint, 2026-11-11 (America/Toronto). Numerics-arm audit, Discovery 3.**

Object audited: Astra, *The completed kernel: a local obstruction and one
exact gluing lemma* (11 November 2026), authoritative copy in the RH goal's
S-matrix lane folder (`Astra_Completed_Kernel_Factor_Map_Obstruction_2026-11-11.md`).
Script: `research/flint_astra_factor_map_audit_2026_11_11.py`; raw output:
`research/flint_astra_factor_map_audit_2026_11_11_output.txt`. Labels are per
item. Nothing here is Lean-checked; the note itself claims paper level only.

**Headline.** The three defects flagged when the note was collected —
a sign/w² mismatch in the prime assembly, a √2 discrepancy in Theorem 5's
constant, an inconsistent archimedean pair — are **not in the authoritative
file**. They were introduced by the on-screen transcription the collector
worked from. In the file, her (7) carries `√(w/2)`, Proposition 2's display
is `−w[q(−a)+q(a)]`, Theorem 5's constant is `√(2w)` with `w = (log 2)/√2`,
and the archimedean maps (6) use `√(m±)` with the pole split (5). Audited
against the file, every flagged point checks out exactly, by hand and by
computation. The obstruction and the CGC reduction stand as paper
mathematics (details and the CGC probe below). CGC itself remains OPEN.

Method: float64 grids; u-grid du = 0.002 (du = 0.0005 in the Theorem-5
window), t-grid dt = 0.02 out to |t| = 100 (dt = 0.05 to |t| = 200 for the
narrow Theorem-5 bump); ψ evaluated by mpmath at 30 dps at grid points.
Refinement pass (du = 0.001, dt = 0.01, |t| ≤ 140) on the headline case
moves ‖P‖² by 2.7e-5 absolute, Q by 1.2e-7 absolute, and the ratio by
9e-8. Quoted ratios are good to about ±2e-7 absolute; Q values to about
±2e-7 absolute for O(1) data. All test functions are real; the bump is
`exp(−1/(1−((u−c)/w)²))` on its support.

## (i) The prime-channel identity — note CORRECT (PROVED by hand; NUMERICAL confirmation)

From her (7), with `a = a_{p,k}`, `w = w_{p,k}`, `f_a(u) = f(u+a)`, on
`I_{R,a} = [−R, R−a]`:

‖P_k f‖² − ‖C_k f‖² = (w/2) ∫ (|f−f_a|² − |f+f_a|²) du
= −2w·Re ∫ f(u)·conj(f(u+a)) du.

Her correlation is `q_{f,g}(u) = ∫ f(v) conj(g(v−u)) dv`, so
`q(−a) = ∫ f(v) conj(f(v+a)) dv` and, substituting `u = v−a`,
`q(a) = conj(∫ f(u) conj(f(u+a)) du)`. Hence `q(a)+q(−a) = 2 Re ∫ f·conj(f_a)`,
and the channel difference is exactly

**−w·[q(a) + q(−a)]**,

which is Proposition 2's display verbatim. No sign error, no w². The w²
suspicion traces to the transcription's `(w/2)[f∓f_a]`: with that prefactor
the difference would be `−(w²/2)[q(a)+q(−a)]`. With the file's `√(w/2)`
prefactor there is no second power of w anywhere.

Numerical confirmation (NUMERICAL). Narrow bump f = bump(0.3, 0.4), R = 2,
channel (p,k) = (2,1): LHS = −0.0000102424; `−w(q+q−)` identical,
relative error 1.9e-13; the sign-flipped candidate `+w(q+q−)` =
+0.0000102424 (killed); the w² candidate `−w²(q+q−)` = −0.0000050201
(off by a factor 2.04, killed). Wide bump f = bump(0, 1), R = 2:

| (p,k) | a | w | LHS | −w(q(a)+q(−a)) | rel. err |
|---|---|---|---|---|---|
| (2,1) | 0.693147 | 0.4901290717 | −0.0688977753 | −0.0688977753 | 0.0 |
| (2,2) | 1.386294 | 0.3465735903 | −0.0037492280 | −0.0037492280 | 1.7e-15 |
| (3,1) | 1.098612 | 0.6342841006 | −0.0307881713 | −0.0307881713 | 1.1e-16 |
| (5,1) | 1.609438 | 0.7197625156 | −0.0006714401 | −0.0006714401 | 6.5e-16 |
| (7,1) | 1.945910 | 0.7354849040 | ~0 | ~0 | both < 1e-10 abs (tail-suppressed overlap sliver) |

Aggregate: for every test function in part (iv), ‖P_R f‖² − ‖C_R f‖² and
Q(f,f) computed independently from her (2) agree to ≤ 2.2e-16 absolute on
the same grid. Proposition 2's signed assembly (8) is an exact identity.

## (ii) Theorem 5's constant — note CORRECT, exact (PROVED by hand; NUMERICAL confirmation)

Window log 2 < 2R < log 3, R = 0.45 in the computation. Here
`w_{2,1} = (log 2)·2^{−1/2} = 0.490129071734`. Take
`I = supp b ⋐ (−R, R−log 2)` — any such I has length < 2R − log 2 < log 3 −
log 2 < log 2 = a, and I, I+a, I−a are pairwise disjoint — and
`f_b = b + b(·−a)`. On `I_{R,a}`: at u ∈ I, f_b(u) = b(u) and
f_b(u+a) = b(u); at every other u ∈ I_{R,a} both values vanish (the second
copy I+a lies above R−a because 2R < log 3 forces R−a < 0 < −R+a). Hence

P_{2,R} f_b = 0,   C_{2,R} f_b = √(w/2)·2b = **√(2w)·b**,
√(2w) = 2^{1/4}·√(log 2) = **0.990079867217**,

which is the file's (14) exactly. The collector's "2wb with w = (log 2)/2"
(coefficient log 2 = 0.6931471806) appears nowhere in the file; the exact
coefficient exceeds it by the factor 1.42838331. There is no √2 discrepancy
in the note to repair.

Numerical confirmation with b = bump(−0.35, 0.08), f_b evaluated exactly
(no interpolation): max|P_{2,R} f_b| = 0.0 and
max|C_{2,R} f_b − √(2w) b| = 0.0 on the grid;
‖b‖² = 0.0106468897 and ‖C_{2,R} f_b‖² = 0.0104367003 = 2w‖b‖² to all
printed digits. Theorem 3's prime-only value on the same witness is
Q_p(f_b, f_b) = −2w‖b‖², with 2w = 0.980258143469 — matching the banked
prime coefficient 0.9802581434685472 from the K=2 reproduction. The
obstruction (P = 0, C ≠ 0, infinite-dimensional negative prime subspace)
stands with the file's constant, unchanged.

## (iii) The archimedean pair (5)/(6) — note CORRECT as written (PROVED by hand; NUMERICAL confirmation)

The file defines, with `r±` the positive/negative parts of a real r,

P_∞ f = (√(m∞⁺)·F, L_c(f)),   C_∞ f = (√(m∞⁻)·F, L_s(f)),
L_c = (L₊+L₋)/√2,   L_s = (L₊−L₋)/√2.

Two lines of algebra: `m∞⁺ − m∞⁻ = m∞`, so the gamma parts of
‖P_∞f‖² − ‖C_∞f‖² give `(1/2π)∫ m∞(t)|F(t)|² dt`; and
`|L_c|² − |L_s|² = L₊·conj(L₋) + L₋·conj(L₊)`, the pole part of her (2).
So ‖P_∞‖² − ‖C_∞‖² is exactly the archimedean-plus-pole summand of Q as
(2) defines it — no swap, no weighted norm, no repair. (Whether one calls
that summand Q_∞ or −Q_∞ is a naming choice for the Weil form's
archimedean piece; the file is internally consistent with the sign in
its (2).)

Numerical confirmation, f = bump(0, 1), R = 2: gamma difference
‖P‖²−‖C‖² restricted to the gamma components = **−0.3058790274**, equal to
`(1/2π)∫ m∞|F|² dt` to all printed digits; pole difference
0.4101108924 = pole part of Q likewise. Side value: m∞(0) = −5.3721834192,
matching the banked archimedean density at 0. The collector's reading —
P and C carrying the multipliers m∞ and |m∞| themselves, not their square
roots — would give `(1/2π)∫ (|m∞|² − m∞²)|F|² dt = 0` identically (m∞ is
real), against the computed −0.3058790274 ≠ 0: that reading is killed,
and it is not the file's. The swapped-roots reading gives +0.3058790274,
wrong sign: also killed.

## (iv) CGC probe (NUMERICAL)

Ratios ‖C_R f‖/‖P_R f‖ under the file's assembly (which (i)–(iii) leave
unchanged). Q(f,f) is R-independent in every row, as (8) requires;
‖P_R f‖² grows with R and the ratio follows her (14a),
ratio² = 1 − Q/‖P_R f‖².

| f | R | ‖P‖² | Q(f,f) | ratio |
|---|---|---|---|---|
| bump(0, 1) | 1.25 | 0.5087762427 | 0.0001252504 | 0.99987690 |
| bump(0, 1) | 2 | 0.7352035689 | 0.0001252504 | 0.99991482 |
| bump(0, 1) | 3 | 1.2107371055 | 0.0001252504 | 0.99994827 |
| bump(0, 1) | 4 | 1.9954160324 | 0.0001252504 | 0.99996861 |
| bump(−0.8, 0.5) | 2 | 0.3302833721 | 0.0009263333 | 0.99859668 |
| bump(−0.8, 0.5) | 3 | 0.5839697052 | 0.0009263333 | 0.99920655 |
| bump(1.2, 0.6) | 2 | 0.5029640853 | 0.0000934226 | 0.99990712 |
| bump(1.2, 0.6) | 3 | 0.8358873192 | 0.0000934226 | 0.99994412 |
| bump(0, 0.35) | 2 | 0.1793738538 | 0.0004185351 | 0.99883266 |
| bump(0, 1.8) | 3 | 2.7014068246 | 0.0000106649 | 0.99999803 |
| autocorr of bump(0, 0.7) | 2 | 0.0353176167 | +1.0e-8 | 0.99999986 |
| autocorr of bump(0, 0.7) | 3 | 0.0584591491 | +1.0e-8 | 0.99999991 |
| autocorr of bump(0, 0.7) | 4 | 0.0966399264 | +1.0e-8 | 0.99999995 |
| f_b, Thm-5 window (unconstrained b) | 0.45 | 0.0276926088 | 0.0090844675 | 0.81972761 |
| f_b, Thm-5 window (moment-constrained b) | 0.45 | 0.0162383379 | 0.0126224948 | 0.47188264 |

Findings:

1. **No ratio exceeds 1 anywhere in the family** — in particular not in
   the f_b directions, the only directions where a single prime channel
   is wholly obstructed. There the completed form is comfortably
   positive: ratio 0.82 for the plain bump witness, 0.47 for the
   moment-constrained witness of Theorem 5 (pole coordinates killed;
   prime-C mass 0.00361583 against gamma-P mass 0.01623834). The
   infinite-rank gamma-to-prime block Theorem 5 forces has, on these
   witnesses, a small load to carry.
2. **The tightest directions are smooth wide data, not f_b.** The
   autocorr-type datum sits at Q = +1.0e-8, i.e. on the CGC boundary
   within this grid's absolute resolution (its cancelling terms are
   ~0.019 in magnitude; ratio 1 − 1.4e-7 against a ratio precision of
   ±2e-7 — its sign is not resolved beyond "≥ 0 within resolution").
   bump(0, 1.8) is next at ratio 1 − 2.0e-6.
3. **Movement with R matches Proposition 6 exactly.** For each fixed f
   the ratio increases monotonically to 1 as R grows, driven entirely by
   ‖P_R f‖² growing at fixed Q — near-cancellation, as she warns, so
   this probe is consistent with CGC and cannot decide it. CGC remains
   OPEN, at exactly ARP strength by her Theorem 4 and her own status
   line.

## Overall verdict

- Proposition 1 (functional-equation folding to (2)): paper proof; read
  in full; the folding, termwise integration on Re s = c > 1, and the
  gamma shift are as stated. Not computer-tested (nothing in it is
  numerical).
- Proposition 2 / (8): exact identity. PROVED by hand above;
  NUMERICALLY confirmed channel by channel and in aggregate.
- Theorem 3 (infinite negative index at one prime, no finite-dimensional
  repair): stands; its witness identity Q_p(f_b,f_c) = −2(log p)p^{−1/2}⟨b,c⟩
  is the (2,1) case priced in (ii).
- Theorem 4 (CGC ⟺ positive factor, with window gluing, dilation and
  shadow actions): read in full; the graph-contraction step, the
  `D_R = (I − T*T)^{1/2}` factor, and the inductive-limit gluing are
  correct as written. The reduction stands.
- Theorem 5 (infinite-rank mixing already at log 2 < 2R < log 3):
  stands with its exact constant √(2w), w = (log 2)/√2.
- Proposition 6 (ratio → 1, no uniform strict contraction): stands;
  confirmed numerically in (iv).
- Controls: F_θ rejected at the local unitary parameter
  (tr V_p = p^θ + p^{−θ} > 2) and Davenport–Heilbronn at n = 6
  (a₆ = 1 ≠ −κ² = a₂a₃) — both are exact algebra in the note, checked
  by reading; neither uses zeros. No zeros are used anywhere in the
  note. Cesàro appears once, as the statehood norm, in its settled form
  sinh(2δT)/(2δT).

**Correction to the collection record** (appended per seat process, not
a correction to Astra): the three flags raised when this note was
collected on 2026-11-11 — prime-assembly sign/w², Theorem-5 √2,
archimedean (5)/(6) inconsistency — came from a defective on-screen
transcription of her (5)–(7) (dropped square roots, a swapped display,
a halved weight). The authoritative file contains none of them.
*Catch it next time:* audit against the authoritative file, never
against a transcription; a "by construction" display that disagrees
with its own definitions is a transcription alarm first and a
mathematics alarm second.

— Flint
