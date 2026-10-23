(* GPP 2026-10-21: prime-specific logarithmic bounded-variation
   bound and global all-time adaptive ultraviolet Fourier tail.
   The analytical theorem uses ψ(x)<=Kx, ψ monotone, Λ(n)>=0;
   no RH or zero archive is used.

   Let q(t)=exp(-t) ψ(exp t)-1, and
   J0(t)=int_0^t exp(u/2)q(u)du.
   On every log block [a,a+1]:
   TV(q)<=3K (positive jumps<=2K, continuous decay<=K).
   h_a(u)=exp(-a/2)J0(a+u) -> h'_a=exp(u/2)q(a+u),
   TV(h'_a)<=sqrt(e)(3K+Q/2)=V, Q=max(1,K-1).
   Subtract the affine endpoint interpolant: g_a(0)=g_a(1)=0.
   Fourier Stieltjes integration twice gives
   |ĝ_a(n)|<=V/(2π² n²) (all nonzero integer n).
   Consequently tail_{|n|>N} <=V²/(6π⁴ N³).
   With N_a=ceil(exp(a/3)(a+1)²),
   Σ_a exp(a)tail_a <=V² π²/5670 <∞.
   This controls UV only; RH-strength low band remains OPEN.
*)
ClearAll[ee,v,n,a,k,kk,x,check,gg];
check[name_,b_]:=If[TrueQ[b],Print["PASS: ",name],Print["FAIL: ",name];Abort[]];
check["zeta-six exactly equals pi^6/945",
 FullSimplify[Zeta[6]==Pi^6/945]];
check["all-time UV constant",
 FullSimplify[(v^2/(6Pi^4))*Zeta[6]==v^2 Pi^2/5670]];
check["Fourier integral comparison",
 FullSimplify[Integrate[x^(-4),{x,k,Infinity}]==1/(3k^3),
 Assumptions->k>0]];
(* Quadratic example with periodic zero endpoints and derivative BV:
   g=u(1-u), TV(g')=2, |g_hat(n)|=1/(2π²n²). *)
gg[u_]:=u(1-u);
check["quadratic Fourier mode bound n=1",
 FullSimplify[Abs[Integrate[gg[x]Exp[-2Pi I x],{x,0,1}]]<=1/Pi^2]];
Print[{"analyticBVlemma",
 "TV(q,[a,a+1])<=3K, from prime jump sum<=2K and decay<=K"}];
Print[{"normalizedPrimeDerivative","d/du (exp(-a/2)J0(a+u))=exp(u/2)q(a+u)"}];
Print[{"alltimeTailBudget","sum_a e^a||P_>Na g_a||²<=V²π²/5670"}];
Print[{"Na","ceil(exp(a/3)(1+a)^2)"}];
Print["No uniform bound of low-band Fourier modes or endpoints follows."];
Print["RH, GRH and zero simplicity OPEN."];
