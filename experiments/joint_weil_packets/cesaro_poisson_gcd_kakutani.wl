(* Golden Physics Project 2026-10-20: exact arithmetic Poisson/GCD
   kernel plus Kakutani critical threshold; NOT an RH proof.
   Product measure is defined coordinatewise for every eps>0.
   For eps<=1/2 it is singular to infinite product Haar, so there
   need not be a global Radon-Nikodym density. *)

ClearAll[theta,r,eps,m,n,check,local,series,aff];
check[label_,x_]:=If[TrueQ[x],Print["PASS: ",label],Print["FAIL: ",label];Abort[]];
local[r_,theta_]:=(1-r^2)/(1-2r Cos[theta]+r^2);
check["every one-prime Poisson density integrates to 1",
 FullSimplify[Integrate[local[r,theta],{theta,0,2Pi}]/(2Pi)==1,
 Assumptions->0<r<1]];
check["critical gcd half-density m12n18",
 FullSimplify[((GCD[12,18]^2)/(12*18))^(1/2)==
 GCD[12,18]/Sqrt[12*18]]];
check["prime-by-prime Fourier product equals gcd kernel",
 FullSimplify[2^(-eps)*3^(-eps)==
 (GCD[12,18]^2/(12*18))^eps,Assumptions->eps>0]];
series=Normal[Series[Sqrt[local[r,theta]],{r,0,4}]];
aff=FullSimplify[Integrate[series,{theta,0,2Pi}]/(2Pi)];
check["local Kakutani Hellinger deficit starts r²/4",
 FullSimplify[Coefficient[aff,r,2]==-1/4]];
Print[{"localHellingerAffinity",aff}];
Print["Kakutani plus prime harmonic divergence gives:"];
Print["eps>1/2 : Poisson product equivalent to Haar"];
Print["0<eps<=1/2 : Poisson product mutually singular with Haar"];
Print["At eps=1/2 covariance EXACT gcd(m,n)/sqrt(mn)."];
Print["At fixed T all prime-power Fourier sums are finite, hence"];
Print["E_eps(T)->E_0(T) as eps->0. But not uniformly in T."];
Print["Elementary PNT gives E_0(T)=o(e^T)."];
Print["For fixed eps>0, mean(P_eps)=o(e^(t/2)) and variance=O(t³),"];
Print["so E_eps(T)/e^T ->4. This is an analytical theorem"];
Print["using prime power independence, NOT verified by this symbolic test."];
Print["Related prior art: Aistleitner/Berkes/Seip (JEMS 2015),"];
Print["GCD sums from Poisson integrals; Kakutani product-measure theorem."];
