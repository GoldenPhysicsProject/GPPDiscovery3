# Per-place and formal-grade census — scope

The BSD-height thinker suggested testing a local/global discrimination
statement with a finite-prime, per-grade signature census. Here we use
the *existing* Wave234 S1 decomposition, not a fabricated local height.

For K=3, prime p, packet b=sin²(πu) supported on [0,1],
a_j=7/20+7j/10, and G_0=b*b^rev, define

 K_{p,theta}(h) = sum_{a>=1} (log p)/p^(a/2) * exp(theta*a*log p)
                                  * G_0(a*log p-h)
 D_ij = K(delta*|i-j|)+K(-delta*|i-j|),
 C_ij = K(2a0+delta*(i+j)+1).
 Q_p^sigma = -(D+sigma*C), sigma = +/-1.

All sums terminate EXACTLY because G0 has compact support.
For shifted F_{1/4}, add the theta=+1/4 and -1/4 matrices.
This isolates only the signed contribution from a *finite prime site* p
inside the full parity form
  P^sigma = Pole^sigma + T_arch - D_prime
                      + sigma*(H_arch - C_prime).
**It is not a full centred form by itself.** In particular,
neither pole nor gamma has been allocated by place. The "grade"
+/-1/2 for Ftheta is FORMAL and fails ordinary motivic integrality.

New numerical check, p=3, K=3, sigma=-1:
  zeta Q_3^- leading minors:
       Delta1≈0.0386849928703
       Delta2≈0.00175354279275
       Delta3≈-0.000149796686154
     inertia (n_+,n_-)=(2,1).
  F_{1/4} combined Q_3^- leading minors:
       Delta1≈0.0857701150863
       Delta2≈0.00850981120310
       Delta3≈+0.000106702337245
     inertia (3,0).
These signs are independently corroborated using mpmath at 90 decimal
digits, but lack directed elementary-function enclosures. This is an
**exploratory sign-flip, not a proof**. An Arb gate should enclose each
entry (exp/log/cos/sin), bound determinants with directed intervals,
then Lean can certify exact rational Sylvester/LDL implications.

**Important conceptual outcome:**
Single prime-power p^a contribution undergoes positive scalar
exp(±theta*a*log p), so its *own* inertia cannot change. But the whole
site sums powers with DIFFERENT positive multipliers, and the site's
matrix inertia CAN differ for finite packets. Therefore claims of
"prime-local blindness" must state whether local means one p^a atom,
one entire prime site's sum, or a local Galois/height complex. None
is a substitute for global zero-localization; the local site
difference does not construct a polarized zeta spectral space.

No BSD or RH theorem is being claimed. If full pole/Arch pairing
changes the classification, that is a distinct calculation.
