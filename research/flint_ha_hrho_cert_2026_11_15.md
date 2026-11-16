# Flint run 16 (2026-11-15): interval certification of (H-A) and (H-ρ)

The last two finite certifications for proof-grade OddFace(0.5), as
catalogued in runs 13/14. Verdicts:

- **(H-A): CERTIFIED WITH A STATED CORRECTION.** The bound holds for
  all j ≥ 73 with certified constant **c_A = 1.7349668350070715**,
  replacing the catalogued float64 value 1.7372355768219303
  (correction −0.13%, filed in §2). The tail side (j > 812) certifies
  at 1.738928806675349, above the catalogued constant; the finite
  range binds, at j = 74.
- **(H-ρ): CERTIFIED AS STATED, with a correction to run 13's filed
  computed value.** ρ_Z ≤ 0.34 is certified; the tightest certified
  threshold on the ladder is **ρ_Z ≤ 0.158** (midpoint float value
  0.15706180861127572). Run 13's filed computed value
  0.33419022381006247 is a float64 artifact and is superseded (§3).

With items 3′ and 4 already banked (run 15; item 4 as GPPVerify
#606), the run-13/14 assembly's named hypotheses are discharged:
(H-Env) by run 14 with Heart's R12 constants, (H-A) and (H-ρ) here.

## 1. Statements and method

L = 0.5, ℓ = 1/4, ω_j = 4πj. St is the bordered Gram after
δ-elimination (basis-independent: the Gram of the δ-projected mode
vectors), X = modes 1–12, Y = 13–72, Z = 73–812,
A = S_Z^eff = St_ZZ − St_ZY St_YY⁻¹ St_YZ.

- (H-A): A_jj ≥ c_A · log ω_j / ω_j⁴ for all j ≥ 73.
- (H-ρ): ρ_Z = ‖D^{-1/2} A D^{-1/2} − I‖ ≤ 0.34, D = diag(A).

Everything is computed in Arb balls (python-flint, prec 200; 𝕂₂
quadrature at prec 160). No float64 result is consumed anywhere;
float64 appears only in midpoint cross-checks, labelled as such.

Pipeline (script §header has the full derivation):
- Modes 1–72, u, g: Heart's Round-7 Arb shards (GPPDiscovery5
  `1da974d`, `experiments/round7_turbo`) — consumed as validation
  targets and as the St_YY source.
- Q(ω) = R(ω) + Π(ω) by Heart's R12 closed form (run 15's code path);
  master 𝕂_j = Q(ω_j)/(πω_j);
  B(ψ_j,ψ_k) = 128(−1)^{j+k}(𝕂_j−𝕂_k)/(ω_jω_k(ω_j²−ω_k²))
  − 2·BI2²·σ_ψj·σ_ψk (j ≠ k).
- Diagonal B(ψ_j,ψ_j) = (128/ω_j²)𝕂₂_j − 2·BI2²·σ_ψj² with
  𝕂₂_j = (1/π)∫_{−ω_j}^{3ω_j} sin⁴(s/2)V(ω_j+s)ℓ²sinc²(sℓ)/(2ω_j+s)² ds
  by panel quadrature (abs_tol 1e-15) plus an explicit tail ball for
  s > 3ω_j. Worst 𝕂₂ relative radius over Z: 0.00114686581218848
  (at j = 74).
- q_δ by run 15's kernel route:
  q_δ = [0.11970390585268734034991818052 ± 6.17e-30], identical to
  run 15's value.
- A assembled by one Arb matrix solve. C_z = zᵀSt_YY⁻¹z (z = q_Y −
  q_δu_Y/g) ≈ 4.0e2 (ball radius 6.85); ‖St_YY⁻¹‖_F ≤
  141268187700.06006; ‖u_Y‖ ≤ 7.266272930744792e-06.

Validations (in-run, all pass; full table in the results JSON):
- V1: the B(ψ,ψ) formula balls lie strictly inside the shard balls
  at (13,14), (13,72), (20,50), (50,72), (31,47) — formula radii
  1e-59…1e-63 against shard radii 2.5e-15…5.4e-14.
- V2: the 𝕂₂ diagonal formula overlaps the shard balls at j = 13, 50,
  72 (formula relative radii 7.9e-3, 1.34e-3, 8.8e-4).

## 2. (H-A)

Finite range (j = 73…812), from diag(A) balls:
profile lower bound min = **1.7349668350070715 at j = 74**
(A_74 = 1.58800935317350324302251190705e-11
± 2.13403717162926974794737167151e-14).
Spot profile lowers: j = 73: 1.84303636290133; j = 812:
2.4181529879570007. Midpoint float profile minimum 1.737301503603918
(run 13 float: 1.7372355768219303).

Tail (all j ≥ 813), analytic lower bound on the profile
A_jj·ω_j⁴/log ω_j from the s-integral form of 𝕂₂ (core |s| ≤ 60 with
the even-envelope parity cancellation of the (1+x)⁻² weight, annulus
60 < |s| ≤ 3000, outer regions by kernel majorants; then
A_jj ≥ B(ψ_j,ψ_j) − ρ̃_j²/g − c_j with c_j bounded through C_z and
‖St_YY⁻¹‖_F; R12 constants C_δ = 19.5771, C_φ = 69.5492):

    profile ≥ 1.738928806675349
    core 1.7701771144837215
    − annulus 0.011039916198792767 − right outer 0.0034981810830465323
    − left1 0.002816480281271141 − left2 0.013882188539570404
    − ρ̃ term 9.71956547509293e-07 − c_j term 1.0569749144205556e-05

Certified: **A_jj ≥ 1.7349668350070715 · log ω_j / ω_j⁴ for all
j ≥ 73** (min of finite 1.7349668350070715 and tail 1.738928806675349).

CORRECTION (H-A). The catalogued c_A = 1.7372355768219303 is run 13's
float64 profile minimum, attained at j = 74, where it equals the
float value to all printed digits — a float minimum cannot certify
at its own value. The interval-certified constant is
1.7349668350070715. Downstream effect (NUMERICAL bookkeeping on run
14's filed components): the tail input (C_y²/c_A)·S_hi moves from
4.900731724582752e-12 to 4.907140201541864e-12 (+0.13%); see §4.

## 3. (H-ρ)

Ã = D^{-1/2}AD^{-1/2} in balls (entry radii: max 0.002688155509531498
on the diagonal, max 0.000278311112197116 off-diagonal). For a
threshold t, ρ_Z ≤ t iff (1+t)I − Ã ⪰ 0 and Ã − (1−t)I ⪰ 0. Each is
certified by a verified LDLᵀ: thin LDLᵀ of the midpoint matrix in
Arb, plus the perturbation bound
‖ |L⁻¹| R |L⁻ᵀ| ‖_∞ < min D with R the entrywise radius matrix. The
factors are dumped (double-double floats, `..._ldlt.npz`); the dump
check re-forms S = M − LDLᵀ in Arb and folds |S| into the bound
(S_max = 1.4760007616132498e-07), so the certificate rests on the
dumped factors, not on the factoring code.

Ladder (certify = dump-check values for t = 0.34 and t = 0.158;
search values otherwise — all agree):

    t       P1=(1+t)I−Ã: dmin / W        P2=Ã−(1−t)I: dmin / W
    0.34    0.3107884620994105 / 0.002708744656308946
            0.3100774538728322 / 0.0027087441149257346     CERTIFIED
    0.25    0.20800477667627948 / 0.002718069405510799
            0.2065745800262493 / 0.002718068289106407      CERTIFIED
    0.20    0.14369091759321448 / 0.002774773848922098
            0.1412327980542723 / 0.0027278061578091288      CERTIFIED
    0.18    0.11431386034794505 / 0.00315081336753555
            0.11100860731473398 / 0.002733606745617771      CERTIFIED
    0.17    0.0981372811551952 / 0.003475923488568016
            0.09415745613080974 / 0.002737149203240044      CERTIFIED
    0.165   0.08856971585736702 / 0.003704799303902817
            0.08507786324890178 / 0.0038369253708326998     CERTIFIED
    0.16    0.07666803659916331 / 0.00400125913575511
            0.07541410214990006 / 0.007265222283092034      CERTIFIED
    0.158   0.07139641733389412 / 0.004949091013751368
            0.06089035996293963 / 0.010945541310639854      CERTIFIED

**ρ_Z ≤ 0.34 CERTIFIED AS STATED; tightest certified ρ_Z ≤ 0.158.**
Midpoint float spectrum of Ã: [0.8429381913887238, 1.1556228568371825],
ρ = 0.15706180861127572 (NUMERICAL cross-check only).

CORRECTION (H-ρ computed value). Run 13 filed ρ_Z =
0.33419022381006247 as a float64 computation. That value is an
artifact, located as follows (all quantities NUMERICAL diagnostics
except the Arb balls):

1. St = M − uuᵀ/g is a cancellation: along run 13's own top
   direction ṽ (D^{-1/2}·top eigenvector of its Ã),
   ṽᵀMṽ = 1138818.4785600835 and (ṽᵀu)²/g = 1138817.1396105096 —
   the St form 1.33894957435697 is the difference of two ~1.14e6
   quantities. A direct Arb evaluation of the same Gram form
   (integral of a square) gives ṽᵀB(φ,φ)ṽ = [1138818.4960 ± 8.04e-5].
2. The routes differ in u: ṽᵀu = −11.873136569080593 (this run,
   formula) vs −11.873132925763487 (run 13 u_all): δ = 3.64e-6,
   worth δ((ṽᵀu)²/g) = 0.366 — the size of the St-form gap
   (1.33894957435697 vs 1.012754108432612 = 0.326, the rest being
   the M-form difference). δSt = −(δu uᵀ + u δuᵀ)/g is a rank-2
   update: it manufactures exactly one spurious eigenvalue.
3. Arb-direct u_j (Heart's R7 evaluator, independent of both
   routes) excludes run 13's u_all at 4 of 8 sampled modes:
   j = 73: deviation −4.492890027058425e-14 vs direct ball radius
   6.32e-15; j = 74: −6.227887989355746e-14 vs 4.70e-15; j = 200:
   +4.695958968441658e-14 vs 2.02e-14; j = 412:
   +3.025020652107302e-14 vs 1.85e-14 (inside at 500, 600, 700,
   812). This run's formula u_j overlapped the direct ball at all
   8. Run 13's u_all came from a kernel-side float formula with
   K′′ interpolation on an 801-point grid (run 13, Part C-lite).
4. Hybrid float test: run 13's St_ZY/St_YY with this run's St_ZZ
   gives ρ = 0.15720037053749525; run 13's St_ZZ gives
   0.33419022381006047. The phantom lives in run 13's St_ZZ, via
   its u_all.

Detector: a spectral quantity of a D-normalised matrix whose
entries are a cancelled difference (M − uuᵀ/g with a 1e6
cancellation along the top direction) must be computed from the
cancelled form in balls; float M and float u separately do not
determine it.

## 4. Downstream (runs 13/14 bounds re-evaluated, NUMERICAL)

Every filed use of ρ_Z was conservative in the safe direction
(1−ρ_Z as a lower-bound factor), so all filed bounds remain valid;
with the certified constants they improve:

- Run 14 tail input: 4.907140201541864e-12 (was 4.900731724582752e-12).
- Assembled v0ᵀLv0 ≤ (1.524894785662644e-11 + 4.907140201541864e-12)
  / (1 − 0.158) = 2.3938346862432667e-11 (was 3.026341802386067e-11);
  room vs β₀ = 1.30163e-8: **543.742643332943×** (was 430.10012913073876×).
- c_dir floor: c_dir ≥ (1−ρ_Z)/RQ_D ≥ 3.5942360966653345e-12
  (was 2.842134835029146e-12 at ρ_Z = 0.33419022381006247);
  718.7200522321574× the threshold c*_v0 = 5.000884677563361e-15.

## 5. Traps caught this run (detectors)

1. **str(arb) is a lossy dump format.** It prints the midpoint to
   ~3 digits and pads the radius to cover the rounding; re-parsing
   inflated 𝕂₂ radii up to ~4× (stored j = 104 radius 5.31e-9 vs
   computed 8.43e-10). Caught by comparing stored balls against
   fresh calls at sample indices. Fix: serialize midpoint and
   radius separately at 30 digits (`ser()` in the script); after
   the fix, stored radii match fresh calls (worst 𝕂₂ relrad
   0.00114686581218848 at j = 74). Detector: after any ball dump,
   round-trip and compare radii at sample indices before trusting
   the dump. (Run 15's committed dump radii are the record and are
   unaffected — they were quoted as committed.)
2. **Block-LDLᵀ transposition slip.** Y = A21·L⁻¹ᵀ requires
   solve(L, A21ᵀ); the first version used solve(Lᵀ, A21ᵀ),
   giving Y = A21·L⁻¹. Pivots alone did not expose it (the
   minimum pivot sat in block 1); the factorisation residual
   S = M − LDLᵀ did (S_max = 0.0082 at entry (184,179)). The
   catalogued 0.34 certificate survived regardless because the
   dump check folds S into the bound; 0.158 needed the fix.
   Detector: always dump-check factors against the matrix
   (residual S); never certify from pivots alone.

## 6. Artifacts and rerun

In GPPDiscovery3, `research/` (branch `flint/moving-in`):
- `flint_ha_hrho_cert_2026_11_15.py` — the certifying script
  (modes: `k2`, `certify`, `verify`).
- `flint_ha_hrho_cert_2026_11_15_k2balls.json` — 𝕂₂ balls, j = 73…812.
- `flint_ha_hrho_cert_2026_11_15_balls.json` — Q balls (j = 1…812),
  A diagonal balls, profile lowers, scalars.
- `flint_ha_hrho_cert_2026_11_15_ldlt.npz` — LDLᵀ factors
  (double-double) for t = 0.34 and t = 0.158, both matrices.
- `flint_ha_hrho_cert_2026_11_15_results.json` — headline numbers;
  this note matches it exactly.
- This note.

Rerun: `python flint_ha_hrho_cert_2026_11_15.py verify` rebuilds A
from the dumped Q/𝕂₂ balls and the R7 shards, re-runs every
certificate and the factor dump check, and compares against the
results JSON. Executed 2026-11-15: all six checks pass
(HA_finite_cA_matches, HA_tail_matches, Hrho_matches,
Hrho_tightest_matches, ldlt_dump_ok, validations_ok).

For the Lean shaping (Claude): (H-A) is a finite list of 740 ball
diagonals plus the analytic tail bound of §2; (H-ρ) is the dumped
L, D factors with the residual+radius perturbation bound of §3 —
the same shape as item 4's export, with L, D as double-double
floats (exact dyadic rationals) in the .npz.

No zero data anywhere. Epistemic status: (H-A) and (H-ρ) are
interval-certified finite computations plus one analytic tail
bound, as labelled per item above.
