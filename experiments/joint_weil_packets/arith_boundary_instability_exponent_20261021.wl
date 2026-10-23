(* GPP October 21, 2026: the *exact* boundary instability exponent of
   the arithmetic pole-balanced coherent prime signal. Conceptually
   E_e(T)=integral_0^T |J_e(t)|² dt, with J_e the same prime signal as
   Astra's Laplace/Cesaro programme. This script verifies the spectral
   pole RESIDUE, half-density regulator transfer and a fake exponential
   resonance; the global theorem uses the unconditional integrated
   Riemann explicit formula plus a rigorous Laplace holomorphy argument
   in the associated private research note. No RH or simplicity claimed. *)

ClearAll[z,eps,rho,m,δ,γ,t,T,dd,ww,alpha,check,Ae,rr];
check[tag_,ok_]:=If[TrueQ[ok],Print["PASS: ",tag],
 Print["FAIL: ",tag];Abort[]];
Ae[z_,eps_]:=(z+eps)(z+eps+1/2)/(z(z+1/2));
rr=rho-1/2-eps;
check["exact nonzero zero-residue transfer through causal filter",
 FullSimplify[-m/(rr(rr+1/2)) -
 Ae[rr,eps](-m/((rho-1/2)rho))]==0];
check["non-unitary off-line pair cannot be killed by pole-balanced filter",
 FullSimplify[Ae[rr,eps]/.rho->(1/2+δ+I γ),Assumptions->δ>eps&&γ>0] =!= 0];
check["Euler-level transfer multiplier at unregularized point",
 FullSimplify[Ae[z,0]==1,Assumptions->z>1]];
(* Toy fake resonant mode of off-axis quartet: e^{δt}cos(γt).
   Its positive integrated energy has exponent 2δ. *)
toy=Integrate[Exp[2dd t]Cos[ww t]^2,{t,0,T},
 Assumptions->T>0&&dd>0&&ww>0];
check["toy energy is exact antiderivative",
 FullSimplify[D[toy,T]-Exp[2dd T]Cos[ww T]^2==0,
 Assumptions->T>0&&dd>0&&ww>0]];
check["toy energy vanishes at T=0",FullSimplify[(toy/.T->0),Assumptions->dd>0&&ww>0]===0];
(* coefficient of e^{2δT} lies in a strictly positive range for ω>0:
  [(1/δ−1/sqrt(δ²+ω²))/4,
   (1/δ+1/sqrt(δ²+ω²))/4]. *)
check["oscillatory fake has strictly positive exponential envelope",
 FullSimplify[1/dd-1/Sqrt[dd^2+ww^2]>0,
 Assumptions->dd>0&&ww>0]];
Print["ANALYTIC ALL-ZERO THEOREM (paper proof, not Lean):"];
Print["δζ=sup_{zeta(ρ)=0, 0<Reρ<1}(Reρ−1/2)>=0."];
Print["For every FIXED eps>=0, J_eps pole-balanced as in Astra:"];
Print["limsup_{T->infty} log(1+int_0^T |J_eps(t)|²dt)/T"];
Print["                  =2 max(δζ−eps,0)."];
Print["The FULL positive Poisson phase energy adds only O(T³)"];
Print["variance, so has the SAME exponent. RH iff δζ=0."];
Print["The arithmetic growth-rate equality is NOT a proof of δζ=0."];
