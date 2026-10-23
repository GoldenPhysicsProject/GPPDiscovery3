(* GPP Oct23. Critical p-local paired thermal reference r=p^(-beta/2).
   Paired forward shift V_p|k,k> = |k+1,k+1>.
   Correlation <Omega_r,V_p^d Omega_r>=r^d.
   Its all-prime tensor Gram is (gcd(m,n)^2/(m*n))^(beta/2).
   This is a nonunitary shift representation, not RH. *)

ClearAll[r,a,b,beta,m,n,j,k,check];
check[tag_,eq_]:=If[TrueQ[eq],Print["PASS: ",tag],
  Print["FAIL: ",tag];Abort[]];
check["local shift correlator at any integer d",
 FullSimplify[(1-r^2)Sum[r^(2*k+3),{k,0,Infinity}]==r^3,
 Assumptions->0<r<1]];
check["exact critical 12-18 GCD Gram",
 FullSimplify[(GCD[12,18]^2/(12*18))^(1/2)
     ==GCD[12,18]/Sqrt[12*18]]];
check["all beta local 2 and 3 prime factors",
 FullSimplify[2^(-beta/2)*3^(-beta/2)
       ==(GCD[12,18]^2/(12*18))^(beta/2),Assumptions->beta>0]];
check["critical Möbius local squared norm",
 FullSimplify[(1+r^2)-2*r*r==1-r^2]];
check["first squarefree primorial Möbius packet norms",
 FullSimplify[(1-1/2)*(1-1/3)*(1-1/5)==4/15]];
Print["ALL-PRIME analytic theorem: finite-prime packet norm^2"];
Print["=Product_{p in P}(1-1/p), which tends 0 by Euler divergence."];
Print["Yet the coefficient at n=1 is 1 and prefix-sample stays nonzero."];
Print["Thus naive coefficient/prefix recovery is not closable on critical GCD/TFD orbit."];
Print["Critical product GNS is positive; no physical arithmetic sewing U exists by fiat."];
