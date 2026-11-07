# From the two-sphere to elliptic arithmetic: exact zeta interfaces

**Scope:** an exact projective-space and elliptic-curve zeta comparison, not a proof of BSD or RH. Suggested by Daniel's question "zeta on S²"; this treats S² as the Riemann sphere P¹(C), with the canonical arithmetic model P¹_Z. This is distinct from the Laplacian spectral zeta on the round sphere and from function composition zeta(s²).

For each prime power q=p^n,
  #P¹(F_q)=q+1.
Thus with T an indeterminate, the point-count definition gives
  Z(P¹/F_p,T)=exp(sum_{n>=1}(p^n+1)T^n/n)
              =1/[(1-T)(1-pT)].
Multiplying local factors (initially where convergent) yields
  zeta(P¹_Z,s)=zeta(s) zeta(s-1).                 (SPHERE)

Let E:y²=x³-x over Q, with p an odd prime of good reduction; a_p=p+1-#E(F_p).
Frobenius acts on H^0,H^1,H^2 with characteristic factors
  P0=1-T; P1=1-a_p T+pT²; P2=1-pT.
Therefore
  Z(E/F_p,T)=P1/(P0 P2).                         (ELLIPTIC)
The H^1 factor is the *new* genus-one channel; it is not Tate-Shafarevich
itself. Its global inverse Euler product defines L(E,s), including separately
specified bad factors. On good primes:
  Z(E/F_p,p^-s) / Z(P¹/F_p,p^-s)
      = 1-a_p p^-s+p^(1-2s)
      = L_p(E,s)^(-1).                            (RATIO)
Hence away from finitely many bad primes, schematically
  zeta_E(s) = zeta(s)zeta(s-1)/L(E,s).
Do not silently convert a good-prime ratio into a global exact identity at
bad fibers; specify a chosen arithmetic model and its correction factors.

For any finite-type arithmetic scheme X and the trivial projective line
bundle P¹_X, the geometric point counts imply the exact identity
  zeta(P¹_X,s)=zeta(X,s)zeta(X,s-1).             (BUNDLE)
In particular for a chosen model X of E, "placing E over S²" produces two
shifted copies of the *same* arithmetic zeta. The two-factor phenomenon is
therefore NOT by itself an off-axis or mixed-state diagnostic.

BSD asserts ord_{s=1}L(E,s)=rank E(Q); its leading term relates a real
period, regulator, Tamagawa factors, torsion and (conditional on finiteness)
#Sha. The algebraic quotient (RATIO) provides a genuine interface between
zeta geometry and the BSD middle factor. It does not prove the rank formula.

Weight-aware mirror: H^0 has Frobenius eigenvalue 1, H^2 has eigenvalue p,
and H^1 has eigenvalues alpha,beta with alpha beta=p. Its Hasse bound
|alpha|=|beta|=sqrt(p) comes from *geometry*, not merely alpha beta=p.
Even valid Hasse bounds on each local Euler factor do not imply the
global analytic continuation's GRH or BSD.

**Exact falsifier: reciprocal duality without positive height.**
Take q=16, J=[[0,1],[-1,0]] and F_bad=diag(8,2),
F_good=[[0,-16],[1,0]]. Both satisfy F^T J F=16 J and
det F=16, hence reciprocal characteristic polynomials.
BAD has eigenvalues 8,2, violating modulus sqrt(16)=4.
For any symmetric H satisfying F_bad^T H F_bad=16H,
H_11=H_22=0; thus det H=-H_12²<=0 and H cannot be
positive definite. GOOD admits H=diag(1,16)>0 with
F_good^T H F_good=16 H. BAD is an ABSTRACT countermodel, NOT Frobenius
of a smooth projective curve: it isolates the added polarization theorem.

If the user means function composition zeta(s²), zeros are preimages
s=±sqrt(rho) of zeros rho of zeta, not a constant-real-part line.
If the user means the round-sphere Laplace spectral zeta, it is
  Z_Delta(w)=sum_{l>=1}(2l+1)[l(l+1)]^-w, Re w>1,
a different object unrelated to the prime Euler product without an
extra arithmetic map.

**Potential next discovery gate:** articulate a global arithmetic pairing
that is (1) built on the actual zeta prime/Arch source, (2) closed in the
critical topology in the presence of the known Möbius nonclosability
counterexamples, (3) polarized/positive on an independently defined
primitive quotient, (4) detects the H^1/defect channel via functorial
comparison with elliptic models while (5) rejecting the Ftheta control.
None of these five load-bearing conditions is proved by (SPHERE)–(BUNDLE).


## Source-only geometric discriminator beyond clock/species counting

For the shifted control F_theta(s)=xi(s-theta)xi(s+theta), ignoring its
Archimedean factors the local Euler product is

  Z_fake,p(T)=1/[(1-p^theta*T)(1-p^(-theta)*T)], T=p^-s.

At p=2, theta=1/4, its first formal point count would have to equal
2^(1/4)+2^(-1/4). Put x=2^(1/4). From 1<x<5/4 (check (5/4)^4>2),
strict AM-GM gives 2<x+1/x, while x+1/x<x+1<9/4<3.
Thus this is NOT an integer, so these local factors cannot be the
ordinary Hasse–Weil point-count zeta of a finite-type variety over F2.
In contrast the two factors of the projective line produce p+1, an
integer for every p. See ftheta_geometric_integrality_gate.py.

This is a *genuine arithmetic/geometric hypothesis* that F_theta fails,
rather than inferring RH from its failure: zeta's own integral prime
factors also do not force zero localization. Stronger future interface:
graded cohomology / effective Frobenius-integrality / pure weights
and a source-derived arithmetic spectral realization, not simply
reciprocal symmetry or 'one clock = one state'. The middle H¹ motive of
E controls its local L-factor; Sha is a distinct global-to-local kernel
in Galois cohomology/ Selmer theory. A proposed BSD-zeta comparison
must construct the map between these structures, not identify them by name.

Relevant prior mathematical framework: Deninger's conjectural
cohomological interpretation of the completed Riemann zeta function
and Nekovar's Selmer/determinant/height complexes. Both are prior art,
not our proof of the requisite zeta cohomology or positive pairing.
