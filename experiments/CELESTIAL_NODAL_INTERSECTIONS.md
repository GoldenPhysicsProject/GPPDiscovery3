# Celestial nodal-intersection theorem: exact scope and falsifiers

Golden Physics Project, 2026-10-25. Extends Daniel Toupin's established October 13 GPP-bridge celestial equator and Hermitian spinor geometry, and the existing GPPVerify ConformalShadowPrincipalSeries.lean. **RH remains OPEN.**

## 1. Dimensional circles coincide rather than intersect discretely

For scalar shadow S_d(Delta)=d−Delta, Re Delta=d/2 is the principal-series axis. Let w_d(Delta)=(Delta−d)/Delta. Then |w_d|=1 iff Re Delta=d/2. Under Delta_d=d*s, w_d(d*s)=(s−1)/s for every s≠0. The d=1,2,3,4 circle images coincide at ALL spectral parameters; unnormalized vertical lines are parallel and never cross. Thus intersection of those circles does not select the discrete zeta zeros. This dimension dictionary is one natural shadow-equivariant CHOICE, not a canonical physical CFT intertwiner.

Choose w=(s−1)/s as stereographic celestial CP1 coordinate, with inverse n(w)=(2Re w,2Im w,1−|w|²)/(1+|w|²). Re s=1/2 corresponds to n_z=0, celestial equator. Anti-linear reflection s→1−conj s sends w→1/conj w and flips n_z. All these facts are already in the project's October 13 notes, not novel work here.

## 2. The EXACT equatorial zero-node intersection theorem

Let xi(s)=s(s−1)π^(−s/2)Γ(s/2)ζ(s)/2 be the classical entire completed xi, G(w)=xi(1/(1−w)) for w≠1. The functional equation xi(s)=xi(1−s), conjugation reality and Cayley reflection imply G(1/conj w)=conj G(w). Hence Im G=0 everywhere on the equator.

**Theorem:** The on-critical-line zeros of xi are precisely equatorial intersections with the real nodal set Re G=0. Every zero of xi, even if off the line, is an intersection of Re G=0 and Im G=0 on the chosen celestial sphere. At a simple equatorial zero the two nodal curves cross orthogonally (holomorphic Cauchy–Riemann), with Im G=0 locally coinciding with the equator. Those other nodal curves are generally NOT literal circles. "All zeros intersect equator" is equivalent to RH; it is not proved by this representation.

**Adversarial fake:** Fθ(s)=xi(s−θ)xi(s+θ), θ=1/4, has the same conjugation and functional equation symmetries but has zeros at rho±1/4 for every known critical-line xi zero rho. Re s=1/4 and3/4 map strictly south/north of the equator. Hence those symmetries and the universal circle do NOT force all zeros to the equator.

## 3. Exact prime-phase interference, not pairwise circle crossing

For 0<Re s<1, eta(s)=sum_(n≥1)(−1)^(n−1)n^(−Re s)exp(−i Im(s)log n)=(1−2^(1−s))ζ(s), so a zeta zero is precisely an infinite cancellation of prime-generated phases exp(−it log n)=product_p exp(−it log p)^v_p(n). This uses the conditionally convergent eta series and its genuine infinite limit, NOT an Euler product convergent in the strip. All prime phases equal 1 at t=0, yet eta(1/2)=0.604898643... !=0. The p=2 and p=3 phases re-align at t=2π/log(3/2)=15.496241677... while the first zeta critical zero is 14.1347251417... . Pairwise phase intersections alone are therefore insufficient. The auxiliary finite-prime phase torus is DISTINCT from Daniel's compact idèle-class K, whose characters grade the arithmetic sectors.

## 4. Exact completed 1D wave integral and replay

Jacobi theta modularity, with H(u)=sum_(n≥1)exp(−π n² exp u), gives for real t:

  xi(1/2+it)=1/2−(t²+1/4) integral_0^infinity H(u)exp(u/4)cos(tu/2)du.

At every critical-line zero the positive-theta-weighted cosine contribution cancels 1/2 exactly. This includes the Archimedean completion and is a classical formula, not a new RH criterion.

Run numerical probe (mpmath optional): python experiments/celestial_nodal_intersections.py --numeric --count 4 --dps 45. First zero t≈14.1347251417346938, w≈0.997500505583138430+0.070659333151277458i (on equator). Integral ≈0.002499494416861569870774151123 = 1/[2(t²+1/4)] to high working precision. Numerically evaluate first four zeros and Fθ off-equator fake zeros. These are numerical checks of established identities only.

## Missing physical/RH mechanisms

Need: zero-independent source-side arithmetic-to-celestial **Hilbert intertwiner**; prime/Arch/pole completion with a nonnegative spectral/source form excluding off-line zeros while rejecting Fθ/DH/Epstein controls; and physical angular/thermal/gravitational state-counting dynamics. Neither shared principal-series circles nor geometric nodal intersections accomplish those. RH and physical 1D unification OPEN.

Prior art: NIST DLMF §25.5.13–14 (theta integral), https://dlmf.nist.gov/25.5 ; Romik's xi Fourier analysis https://math.ucdavis.edu/~romik/data/uploads/papers/riemannxi.pdf ; project October 13 celestial notes.
