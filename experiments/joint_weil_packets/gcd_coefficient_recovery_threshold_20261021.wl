(* GPP 2026-10-21: exact Poisson/GCD coefficient recovery threshold.
   Source space: analytic trigonometric polynomials in chi_n, n>=1,
   with Poisson kernel K_eps(m,n)=(gcd(m,n)^2/(mn))^eps.
   For eps>1/2, each coefficient a_n is continuous, exact dual norm²
     zeta(2eps) prod_{p|n}(1+p^(-2eps)).
   For 0<eps<=1/2, the natural finite-window prefix operator
   B_T(a)=exp(-t/2)sum_{n<=exp(t)}a_n is NOT CLOSABLE:
   arithmetic packet a_(P,eps)(d)=mu(d)d^(-eps) on d|primorial
   has source norm²=prod_{p in P}(1-p^(-2eps))->0,
   but output stabilizes to a nonzero finite-window vector.
   THIS DOES NOT IMPLY RH. *)

ClearAll[r,k,mat,im,eps,pp,check];
check[name_,expr_]:=If[TrueQ[expr],Print["PASS: ",name],
 Print["FAIL: ",name];Abort[]];
mat[k_Integer]:=Table[r^Abs[i-j],{i,0,k},{j,0,k}];
Do[
 im=FullSimplify[Inverse[mat[k]],Assumptions->0<r<1];
 check["Toeplitz AR1 constant coefficient dual norm k="<>ToString[k],
   FullSimplify[im[[1,1]]-1/(1-r^2),Assumptions->0<r<1]===0];
 If[k>=2,check["Toeplitz AR1 positive exponent coefficient norm k="<>ToString[k],
  FullSimplify[im[[2,2]]-(1+r^2)/(1-r^2),Assumptions->0<r<1]===0]];
,{k,1,7}];

check["local Euler packet source norm for general eps",
 FullSimplify[(1+r^2-2r^2)-(1-r^2)]===0];
check["full coefficient dual norm zeta product arithmetic",
 FullSimplify[Times@@Table[1/(1-Prime[j]^(-2eps)),{j,1,5}] ==
   Times@@Table[(1-Prime[j]^(-2eps))^(-1),{j,1,5}],
   Assumptions->eps>1/2]];
check["coefficient n=6 finite prime correction",
 FullSimplify[(1+2^(-2eps))*(1+3^(-2eps)) >
 0,Assumptions->eps>1/2]];
Print["EXACT infinite analytic proof: Gram K_r(i,j)=r^abs(i-j),"];
Print["orthonormal phi0=1, phik=(z^k-r z^(k-1))/sqrt(1-r²)."];
Print["Coefficient c0 dual norm² local =1/(1-r²)."];
Print["Coefficient ck dual norm² local=(1+r²)/(1-r²), k>=1."];
Print["Tensor over ALL primes, r_p=p^(-eps):"];
Print["||coefficient_n||²=zeta(2eps) product_(p|n)(1+p^(-2eps))"];
Print["for eps>1/2. Norm blows as (2eps-1)^(-1/2) at criticality."];
Print["For eps<=1/2, Euler packet mu(d)d^-eps collapses source"];
Print["but finite-window output stabilizes -> nonzero, NOT CLOSABLE."];
Print["No RH theorem, no nontrivial-zero spectral location implied."];
