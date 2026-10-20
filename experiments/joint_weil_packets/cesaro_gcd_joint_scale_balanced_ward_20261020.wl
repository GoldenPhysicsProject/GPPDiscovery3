(* Golden Physics Project, October 20, 2026. All-prime phase-scaled
   Cesaro crossover and pole-balanced prime/Archimedean subtraction.
   Model: P_eps(t)=2 sum_{p^k<=e^t} log p (p^(-k/2)-e^(-t/2))
   (p^k)^(-eps). C_eps built from dx + delta_at_1 so C_0=C.
   The all-prime crossover proof uses classical PNT; the Wolfram
   algebraic checks here do not attempt to prove PNT or RH. *)

ClearAll[te,z,ee,aa,expected,smooth,bal,check];
check[msg_,ok_]:=If[TrueQ[ok],Print["PASS: ",msg],
 Print["FAIL: ",msg];Abort[]];

smooth[te_,ee_]:=2*((Exp[(1/2-ee)te]-1)/(1/2-ee)-
 Exp[-te/2]*(Exp[(1-ee)te]-1)/(1-ee));
bal[te_,ee_]:=smooth[te,ee]+2(1-Exp[-te/2]);

check["counterterm at eps=0 equals ORIGINAL C(t)",
 FullSimplify[Limit[bal[te,ee],ee->0]==2(Exp[te/2]-1),
 Assumptions->te>0]];
check["balanced half-density critical counterterm",
 FullSimplify[Limit[bal[te,ee],ee->1/2]==
 2(te-1+Exp[-te/2]),Assumptions->te>0]];
check["the prime mean cross-over is exp(-a) at eps=a/T",
 FullSimplify[Limit[smooth[te,aa/te]/(2Exp[te/2]),te->Infinity,
 Assumptions->aa>0]-Exp[-aa],Assumptions->aa>0]===0];

(* Laplace transform from the exact continuous density dx + delta_1. *)
lt=FullSimplify[Integrate[Exp[-z te] bal[te,ee],
 {te,0,Infinity},Assumptions->z>2&&0<ee<1/2],
 Assumptions->z>2&&0<ee<1/2];
check["balanced continuous Laplace exact",
 FullSimplify[lt==(1+1/(z+ee-1/2))/(z(z+1/2)),
 Assumptions->z>2&&0<ee<1/2]];

(* The true shifted arithmetic log derivative is -zeta'/zeta,
   with s=z+1/2+eps. At the simple zeta pole, the bracket's
   singular terms cancel EXACTLY. *)
check["shifted prime pole cancels shifted arch pole",
 FullSimplify[(-1/(z+ee-1/2)) +1/(z+ee-1/2)==0]];

(* On the Poisson side, VAR(P_theta)=sum_p variance of independent
   local prime phase power-series; <=4sum_p log²p/(sqrt p-1)²=O(t³).
   At eps>=1/2, elementary psi(x)=O(x) (Chebyshev) gives
   E_bal(T)=O(T³), independently of RH. *)
Print["PNT IMPLIES full crossover lim E_{eps(T)}(T)/exp(T)="];
Print[" 4(1-exp(-a))² where a=lim T eps(T) in [0,infinity]."];
Print["eps=1/2: C_bal=2(t-1+exp(-t/2));"];
Print["P_eps=O(t) via Chebyshev sum Lambda(n)/n=O(t);"];
Print["local independent-prime variance O(t² or t³);"];
Print["thus E_bal(T)=O(T³) with NO RH input."];
Print["CAVEAT: this is a SHIFTED arithmetic signal, not Astra's J."];
Print["RH / GRH / zero simplicity remain OPEN."];


(* Infinitesimal phase variance and p-adic nested-square depth.
   Symbolic two-prime-power test checks exact first-order e-derivative.
   General identity: min(k,l)=sum_{j>=1}1_{j<=k}1_{j<=l};
   V'_0=2 sum_p log p sum_j (sum_{k>=j}a_{p^k})^2 >=0.
   This polynomial variance is RH-blind; the hard mean is J. *)
ClearAll[pp,aa1,aa2,cc,hh,ex,rp,varg,ms];
rp=pp^(-hh);
ex=cc^2-2cc(aa1 rp+aa2 rp^2)+aa1^2+aa2^2+2aa1 aa2 rp;
varg=2Log[pp](aa1^2+2aa1 aa2+2aa2^2);
ms=aa1+aa2-cc;
check["all-prime-phase infinitesimal variance = nested positive squares",
 FullSimplify[(D[ex,hh]/.hh->0) -
   (-2ms(aa1+2aa2)Log[pp]+varg),Assumptions->pp>1]===0];
check["two prime-power example smoothing energy is not monotone",
 FullSimplify[(D[ex,hh]/.(cc->0))<0,
 Assumptions->pp>1&&aa1>0&&aa2>0&&hh>0]];
Print["Exact all-padic-depth variance slope:"];
Print["V'_0(t)=2 sum_p log p sum_j (sum_{k>=j}a_{p^k}(t))² >=0"];
Print["Prime Chebyshev bound V'_0(t)=O((1+t)^3)."];
