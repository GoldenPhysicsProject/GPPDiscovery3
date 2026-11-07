# Massive-vector three-polarization sewing: corrected 5D Ward reduction

Codex/GPT discovery track, 2026-08-26.

## Correction to the earlier shortcut

The earlier version of this note proposed reducing the massive projector by proving the four-dimensional transversality conditions

\[
p\cdot J_{4D}=0.
\]

That condition is **not generically true** for the four-dimensional massive vector obtained from a five-dimensional massless gauge field.  The executable audit in
`massive_vector_5d_ward_reconstruction.py` verifies the full five-dimensional Ward identity while finding nonzero four-dimensional longitudinal contractions.

Write the null five-dimensional cut momentum as

\[
K^A=(p^\mu,\kappa),\qquad K^2=0,\qquad p^2=\kappa^2=\mu^2.
\]

For a five-dimensional gauge current \(J^A=(J^\mu,J^5)\), when all other exposed vector legs are in physical states, the correct Ward identity is

\[
K_AJ^A=0
\quad\Longrightarrow\quad
p_\mu J^\mu=\kappa J^5
\]

in mostly-minus signature.  Thus the longitudinal part of the four-dimensional massive-vector completeness projector does **not** vanish.  Instead,

\[
\frac{(p\cdot J_L)(p\cdot J_R)}{\mu^2}
=J_L^5J_R^5
\]

when the same signed fifth momentum is used in the sewn orientation. Therefore for a single exposed physical vector leg,

\[
\boxed{
-J_{L,4}\cdot J_{R,4}
+\frac{(p\cdot J_{L,4})(p\cdot J_{R,4})}{\mu^2}
=
-J_{L,4}\cdot J_{R,4}+J_L^5J_R^5
=
-J_{L,5D}\cdot J_{R,5D}.
}
\]

No high-energy or small-`mu` approximation is involved.

## Exact projector identity

For a four-dimensional massive cut momentum `p` with `p^2 = mu^2 > 0`, the sum over the three physical vector polarizations is

\[
\sum_{\lambda=1}^{3}\varepsilon^{(\lambda)}_\mu(p)\,
\varepsilon^{(\lambda)*}_\nu(p)
=
-\eta_{\mu\nu}+\frac{p_\mu p_\nu}{\mu^2}
\]

in mostly-minus conventions.

Hence for two exposed four-dimensional currents,

\[
\mathcal S_V
=
-J_{L,4}\!\cdot J_{R,4}
+\frac{(p\!\cdot J_{L,4})(p\!\cdot J_{R,4})}{\mu^2}.
\]

For one exposed vector leg with every other vector leg physical, the five-dimensional Ward identity converts this exactly to

\[
\boxed{\mathcal S_V=-J_{L,5D}\cdot J_{R,5D}.}
\]

## Second correction: do not double-metric-contract the raw rank-two tensor

The two-particle cut exposes two massive-vector legs simultaneously.  A tempting further shortcut is to form an unprojected rank-two tree tensor \(J^{AB}\), assert

\[
K_{1A}J^{AB}=0,
\qquad
K_{2B}J^{AB}=0,
\]

for arbitrary values of the other exposed index, and replace both physical projectors by bare five-dimensional metrics.  The new executable audit
`double_massive_vector_projector_audit.py` shows that this stronger tensor statement is false in general: when the second exposed leg is filled by an arbitrary basis vector rather than a physical transverse polarization, the first Ward contraction is generically nonzero, and vice versa.

This does **not** violate gauge invariance.  The ordinary on-shell Ward identity assumes the other external vector states are physical.  It means only that the Ward reconstruction may not be applied twice naively to the completely unprojected rank-two tensor.

The safe exact two-vector object is therefore the physical double projector

\[
\boxed{
J_L^{\mu\nu}
P^{(1)}_{\mu\rho}
P^{(2)}_{\nu\sigma}
J_R^{\rho\sigma},
\qquad
P^{(i)}_{\mu\nu}
=-\eta_{\mu\nu}+\frac{p_{i\mu}p_{i\nu}}{\mu_i^2}.
}
\]

The new audit independently evaluates the explicit \(3\times3\) polarization sum and this double-projector contraction and finds equality to floating-point precision over deterministic kinematic trials.  Thus the massive projector is not merely formal bookkeeping: it is the honest compact representation of the nine physical vector states.

## Consequence for the `D_s=4`, nonzero-`mu` baseline

The previously established reconstruction remains

\[
C^{(4)}(\mu)=C^{(V_m)}(\mu)-C^{(S)}(\mu).
\]

For `C^(V_m)` with two internal vector lines, retain both massive projectors explicitly (or perform the equivalent explicit physical-polarization sum).  The fifth-current reconstruction remains useful for a single leg after the other vector states are physical, but a raw double five-dimensional metric contraction is not justified solely from the unprojected rank-two tensor.

## What remains for the MHV cut

The decisive calculation is now:

1. construct the complete color-ordered five-dimensional/KK four-gluon tree with the two massive cut momenta;
2. retain the two honest four-dimensional massive projectors while sewing the two trees, or equivalently sum the three physical states on each line;
3. simplify that projected tensor expression using Ward identities only after the companion exposed vector state is physical/projected;
4. subtract `C^(S)` to obtain the honest nonzero-`mu` `D_s=4` baseline;
5. restore general `D_s` through
   \[
   C^{(D_s)}=C^{(V_m)}+(D_s-5)C^{(S)}.
   \]

Only after this one-loop numerator is honest should the same projector-level state-sum discipline be exported to gravity and higher/generalized cuts.

## Verified discovery evidence

`massive_vector_5d_ward_reconstruction.py` builds the complete color-ordered five-dimensional four-gluon tree from the cubic and quartic Yang-Mills vertices. On deterministic generic transverse states it verifies the five-dimensional Ward replacements while explicitly finding nonzero `p.J_4D`, and it verifies

\[
p\cdot J_{4D}-\kappa J^5=0.
\]

`double_massive_vector_projector_audit.py` then exposes both vector legs.  It records the negative result that the raw rank-two tensor is not separately transverse against an arbitrary unphysical basis choice on the companion leg, and verifies the positive result

\[
\text{double massive-projector contraction}
=
\text{explicit nine-state physical-polarization sum}.
\]

These are reproducible discovery results and exact algebraic guides.  The final MHV cut numerator is not yet promoted to Lean.
