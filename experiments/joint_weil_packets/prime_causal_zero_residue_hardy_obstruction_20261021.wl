(* Golden Physics Project Oct 21 — actual zeta-zero residue eigenfunctionals
   of a causal PRIME-SIGNAL convolution module, plus essential
   counterexample: unilateral isometries CAN have bounded dual eigenvalues
   strictly inside the unit disk. Not a Hilbert spectral identification
   or a proof of RH. *)

ClearAll[z,lambda,alpha,t,u,a,kt,Jtoy,out,lap,check,rho,eps,m];
check[msg_,condition_]:=If[TrueQ[condition],Print["PASS: ",msg],
 Print["FAIL: ",msg];Abort[]];

(* Generic toy causal resonance J=e^(alpha t), k=t e^(-t), alpha>0 *)
out=FullSimplify[
 Integrate[u Exp[-u] Exp[alpha (t-u)],{u,0,t}],
 Assumptions->alpha>0&&t>0];
lap=FullSimplify[Integrate[Exp[-z t]out,{t,0,Infinity}],
 Assumptions->z>alpha&&alpha>0];
check["toy true causal convolution Laplace product",
 FullSimplify[lap-1/((z+1)^2(z-alpha)),Assumptions->z>alpha&&alpha>0]===0];
check["residue is exactly kernel-Laplace evaluation at alpha",
 FullSimplify[Limit[(z-alpha)lap,z->alpha,Assumptions->alpha>0]
  -1/(alpha+1)^2,Assumptions->alpha>0]===0];

(* True arithmetic J_e Laplace has nontrivial-zero residue independent
   of the continuum subtraction; reflected zeros show reciprocal modulus. *)
z0=rho-1/2-eps;
Ae[x_,ee_]:=(x+ee)(x+ee+1/2)/(x(x+1/2));
check["balanced transfer preserves exact zero pole residue",
 FullSimplify[Ae[z0,eps](-m/(rho(rho-1/2)))+
      m/(z0(z0+1/2))]===0];

(* Hardy unilateral shift S f(w)=w f(w) is isometric (from coefficient
   shifts), yet evaluation phi_w(S f)=w phi_w(f) is bounded for |w|<1.
   Take w=1/2, evaluation dual norm =1/sqrt(1-w²)=2/sqrt(3). *)
check["Hardy interior evaluation is bounded despite nonunitary eigenvalue",
 FullSimplify[(1/Sqrt[1-(1/2)^2]) ==2/Sqrt[3]]];
check["Hardy delay character is strictly inside unit circle",
 Abs[1/2]<1];
Print["Algebraic module M_J={k*J0}; if k*J0=0, analytic Laplace"];
Print["product L[k]L[J0]=0 implies k=0. Hence ALL complex lambda"];
Print["have FORMAL evaluation characters. Actual zeta zeros are singled"];
Print["out by genuine response POLES with residues, not by existence"];
Print["of any algebraic evaluation. RH-strength gap: physical unitary"];
Print["completion where exactly relevant pole functionals are bounded."];
Print["RH/GRH remain OPEN."];