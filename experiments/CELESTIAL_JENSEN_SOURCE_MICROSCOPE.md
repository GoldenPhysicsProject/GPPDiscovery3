# Celestial local Jensen microscope: a zero-independent completed-theta source and rigorous control separation

Golden Physics Project — Codex discovery, 2026-10-25, continuing Daniel Toupin's Holographic Chain, the preexisting Oct 13 celestial Cayley geometry and merged Discovery2 PR24. **RH NOT PROVED.** This is a local source-side diagnostic / adversarial falsifier; do NOT present Jensen's theorem as progress on the missing positivity estimate.

## 1. Exact local celestial disk charts, not merely a plotted intersection

For h in R and a>0 set
   z_(h,a)(s)=(s−1/2−ih−a)/(s−1/2−ih+a);
   s_(h,a)(z)=1/2+ih+a(1+z)/(1−z).
Then |z|<1 iff Re(s)>1/2; |z|=1 iff Re(s)=1/2 (away from the chart pole). Each chart identifies the RIGHT half-plane with one celestial hemisphere, and its boundary with the SAME celestial equator; changing h,a is a boundary-preserving Möbius reparametrization, NOT a new physical celestial sphere or a dynamical dimensional intertwiner. At s=1/2+ih+δ with δ>0, z=(δ−a)/(δ+a). Compared with default (s−1)/s, choosing h near the height of interest moves very-high zeros far from the near-singular stereographic point z=1; it supplies a local spectral microscope.

The SO(3)/round-sphere geometric latitude is n_z=(1−|z|²)/(1+|z|²). It is positive for the right hemisphere. Möbius changes of chart generally are conformal maps of the sphere, NOT SO(3) rotations; the physical interpretation of these points requires an actual CFT action/intertwiner.

## 2. Jensen and argument-principle theorem, exact scope

Let A(s) be analytic throughout Re(s)>1/2 with entire continuation as needed, and A(s_(h,a)(0))≠0. Set g(z)=A(s_(h,a)(z)) on |z|<1. For any 0<r<1 with no zero on |z|=r, Jensen's classical theorem gives the NONNEGATIVE defect
  J_A(h,a;r):=(1/2π)∫_0^(2π) log|g(re^(iθ))/g(0)| dθ
             = Σ_(|z_k|<r) m_k log(r/|z_k|) ≥0.
Zero at chart center would require factoring z^m; choose nonvanishing center. The argument principle gives wind_(|z|=r) g=Σ_(|z_k|<r)m_k.

For A(s)=xi(s), the completed Riemann xi has its zeros only in 0<Re(s)<1, and the functional equation reflects zeros across Re(s)=1/2. Therefore RH is equivalent to J_xi(h,a;r)=0 for **every** center (one h,a already covers the entire right half-plane as r↑1). This is just an equivalence; the proof barrier is to derive the vanishing Jensen defect from independent prime–Arch source estimates, not to assume it. A finite numerical winding 0 cannot certify the entire infinite half-plane.

## 3. Exact shifted Fθ discrimination at a near-critical height

Let Fθ(s)=xi(s−θ)xi(s+θ), θ=1/4. Take a=1/2 and h=γ where ρ=1/2+iγ is any (known to exist) critical xi zero. Then s*=3/4+iγ is a zero of Fθ and is mapped to
   z_(γ,1/2)(s*)=(θ−a)/(θ+a)=−1/3.
This is an **exact** location conditional only on the already known existence of ρ (Hardy), with NO numerical search. At r=9/20=0.45, its guaranteed Jensen contribution is
   J_Fθ(γ,1/2;9/20) ≥ log((9/20)/(1/3))=log(27/20)
                                        =0.300104592450338...
One included zero suffices; any additional zeros only increase this number. No cancellation possible because each Jensen term is nonnegative. Winding is at least one, and the exact detection radius r>1/3 is independent of γ after recentering.

In the numerical finite contour centered at mpmath's first critical zero γ≈14.13472514173469 (NOT independently interval-certified by this script), samples=192, dps=32 gave xi winding 0 and fake winding 1 at r=0.45. Numerically J_xi≈0 and J_Fθ≈log(27/20), consistent with exactly one zero of fake. These numerical contour estimates are not proof of the respective counts absent certified nonvanishing arcs and phase bounds.

Notice the fake is even nonnegative **on the actual equator**: for real t, by xi functional equation and reality,
   Fθ(1/2+it)=|xi(1/2+θ+it)|² ≥0.
A pointwise positive equatorial wave can carry off-equator zeros! In contrast real xi(.5+it) changes sign at simple zeros.

## 4. Entire completed theta source: no continuation trick inside the contour

Let H(u)=Σ_(n≥1)exp(−πn² e^u), u≥0. Classical Jacobi transformation gives the following entire identity for every complex s:
  xi(s)=1/2+s(s−1)∫_0^∞ H(u)e^(u/4)cosh((s−1/2)u/2)du.
For s=.5+it this becomes the previously banked real cosine interference formula. The super-exponential tail provides absolute convergence uniformly on compact subsets of C, so contour source values are well defined without any zero list. Numerical mpmath evaluations of the source integral at s=.75+iγ and s=1+iγ match direct gamma-zeta xi at ≲1e−32 absolute with 32 digits, but this check is numerical only. SOURCE includes full Archimedean completion and the theta lattice; it is not individually separated prime factors.

## 5. Positive genuine theta Fourier density STILL does not force the equator

A known integral-by-parts refinement (Riemann's Fourier kernel; also prior team θ-control) gives
  xi(1/2+it)=∫_R Φ(x)e^(itx)dx, where Φ(x)=Φ(−x) and for x≥0
  Φ(x)=Σ_(n≥1)[4π² n⁴ e^(9x/2)−6πn² e^(5x/2)]e^(−πn²e^(2x))
      =Σ_(n≥1)2πn²e^(5x/2)(2πn²e^(2x)−3)e^(−πn²e^(2x)) >0.
The positivity is termwise, as 2πn²e^(2x)>3. Exponential tilts are integrable at every finite real θ thanks to the superexponential theta tails. For shifted comparison define
  Ψθ(u)=∫_R Φ(x)Φ(u−x)cosh(θ(u−2x))dx >0.
It is even and decays enough to define an entire Fourier transform, and
  Fθ(1/2+it)=∫_R Ψθ(u)e^(itu)du.
This is exactly the previous programme's tilted convolution Φθ=(Φe^(−θx))*(Φe^(+θx)), rewritten symmetrically. **Prior art/project status:** team's October 1 θ-control already had the convolution, not a new discovery. The new observation here combines its strict positive even entire density and Bochner positive-definiteness WITH a localized celestial Jensen count that detects the explicit off-equator zero.

So even (i) positive theta-wave density; (ii) positive-definite real-axis amplitude; (iii) full Archimedean completion; (iv) modular theta functional equation; (v) an exact celestial anti-linear equator mirror; (vi) pointwise nonnegative equatorial wave DO NOT imply all nodes are equatorial. A genuinely arithmetic, unshifted Euler-product/weight/physical-sector discriminator or full Weil inequality is still required, independently verified rather than renamed.

## 6. Tests and proof gap

Code experiments/celestial_jensen_source_microscope.py; tests/test_celestial_jensen_source_microscope.py. Stdlib exactly checks Cayley algebra, polynomial zero winding/Jensen on a sample circle, and exact fake root −1/3. With optional mpmath --numeric --dps 38 --samples 160 replays known gamma, xi and Fθ contour integrals, two independent completed source theta anchors, Φ positivity and Ψθ convolution samples. Tests are not arbitrary interval certification; numerical zero counts are exploratory only. All exact theorems in the note follow from classical Fourier/Jensen/functional equation and explicit algebra, not from the finite numerical samples.

**Next actual obstruction:** Find a *zero-independent*, quantitative prime–Archimedean control ensuring log|xi(s_(h,a)(z))| is harmonic for every disk in the right half-plane (J=0), while the same method fails strictly for Fθ, Davenport–Heilbronn, and Epstein forms. That statement is equivalent in conclusion to RH, not a solution by itself; its *source-side proof* is missing. The local chart makes the failure visible at O(1) disk radius instead of radial resolution ~γ^(-2), allowing fair adversarial numerical/certified testing.
