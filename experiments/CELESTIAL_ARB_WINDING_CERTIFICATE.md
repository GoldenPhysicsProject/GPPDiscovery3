# Certified local celestial winding from full-arc Arb enclosures

Golden Physics Project, 2026-10-25. Continues the established Oct13 celestial Cayley geometry and Discovery2 PR24/PR25. **Finite-local certified statement, not RH.**

## The exact object

Use the height-adapted Cayley inverse
  s(w)=1/2+ih+(1/2)(1+w)/(1−w),
where the height h is the *exact decimal rational* 14.1347251417346937904572519835624702707842571156992431756855. Arb encloses this rational, not a guessed floating-point center. Take the circle |w|=r=9/20. Since r<1, the whole closed disk stays in Re s>1/2. Its explicit Re s range is
  1/2 + (1/2)(1−r)/(1+r) <= Re s <= 1/2 + (1/2)(1+r)/(1−r).
At r=9/20 this is exactly
  1/2+11/58=20/29 <= Re s <= 1/2+29/22=20/11.
The imaginary coordinate is centered near the first known critical zero height. The completed ENTIRE function xi(s)=s(s−1)pi^(−s/2)Gamma(s/2)zeta(s)/2 is evaluated via Arb interval special functions (not zero-list interpolation). Control F_(1/4)(s)=xi(s−1/4)xi(s+1/4).

The first critical zero rho is separately enclosed using python-flint Arb's zeta_zeros(1,1), which has existing Grok S51 independent use; the height h differs by <1e−30 from the certified gamma enclosure. Its shifted control zero lies near local w=−1/3, but this proximity is not used to count zeros: we independently interval-evaluate the COMPLETE analytic contour.

## Why the count is rigorous, unlike the earlier numerical PR25

Partition theta∈[0,2π] into N=2048 exact adjacent intervals I_k=[2πk/N,2π(k+1)/N]. Build theta_k as an Arb interval with midpoint π(2k+1)/N and radius π/N. It encloses the ENTIRE segment, not just one sample. Evaluate w=r exp(i theta_k), then s(w), then both xi(s(w)) and Ftheta(s(w)) with directed acb rectangular complex enclosures.

**For every k and both functions**, the interval image excludes 0. The smallest certified lower absolute value across all full arcs was
- xi: 0.00025157327598209985140 (strictly positive directed lower bound);
- Ftheta: 3.5246467504064057861e−8 (strictly positive directed lower bound).
These are interval lower bounds, NOT apparent sample minima. Consequently each true complex image arc lies within a zero-free convex rectangular acb enclosure. Any such image arc is homotopic, with endpoints fixed, to the straight segment joining its endpoint images inside that convex rectangle and never crossing the origin. Replacing all arcs by their chords therefore preserves the exact integer winding number.

Evaluate the complex-ball endpoint ratios f(s(w_{k+1}))/f(s(w_k)), take Arb interval arguments, and prove every phase increment lies strictly in (−π/2,+π/2). Summing the directed intervals and dividing by 2π yields:
  xi turns ∈ [+/- 2.64e−70], uniquely ENFORCING integer winding 0;
  Ftheta turns ∈ [1.00000000000000000000000000000 +/- 3e−35], uniquely ENFORCING integer winding 1.
The upper bound on individual phase increments was ~0.00474795 radians (xi), ~0.00976372 radians (Ftheta). These small increments are certified intervals. Analytic argument principle now gives:

**Theorem (finite celestial disk zero count, independently Arb-certified):** xi(s) has exactly ZERO zeros (counting multiplicity) in the specified disk |w_(h,1/2)(s)|<9/20, while F_(1/4)(s) has exactly ONE zero (counting multiplicity) there.

This is not merely an RH-equivalent rewrite: one actual explicitly described disk is certified and the control's off-equator zero is rigorously detected. It does not put any bound on zeros outside this disk or at another height. In particular, it does NOT approach a global RH proof without an independent all-heights source inequality.

## Programmatic evidence

Script experiments/celestial_arb_full_arc_winding.py, reproducible with python-flint==0.9.0, segments=2048, bits=256. CI Arb crosscheck printed RIGOROUS_ARB_COMPLETE_ARC_WINDING_CERTIFICATE at workflow run 40166044252, job 121087190741, and succeeded. Three Python stdlib jobs also green at time of this note; source lane should be awaited before merging. The script aborts upon ANY undecidable/nonzero-overlapping arc ball, failed gamma enclosure, inconclusive phase interval or non-unique winding integer. No extrapolation from finite sampled points is allowed.

The certified region's quoted σ bound is just the image of this particular disk; it is NOT a general zero-free strip. This is an interval result about analytic xi, NOT a mathematical proof from first-principles prime source of the full Riemann hypothesis.

## Why this matters for Daniel's programme

- Both functions inherit an equator and an exact Hermitian/reflection symmetry, and both have a completed theta wave; the shifted control even has a strictly positive Fourier convolution source, per the team's earlier Oct1 theta-control.
- Nevertheless the local celestial zero-vorticity counts DIFFER. This demonstrates that a valid source discriminator must encode more than representation geometry, Bochner positivity, spectral principal lines and positive boundary amplitudes.
- The height-centered chart removes an artificial numerical resolution obstacle. It does not supply physical dynamical intertwiners between d-dimensional CFTs; celestial CP1 here is a chosen stereographic sphere.
- For any proposed prime–Arch weighted Hilbert/polarization construction, the operational acceptance gate is: explain from genuine arithmetic data why every such right-hemisphere disk has zero winding for xi, while an analogous attempted proof FAILS at the certified Ftheta disk and other non-Euler adversaries. Such a universal proof is still missing. RH/BSD OPEN.
