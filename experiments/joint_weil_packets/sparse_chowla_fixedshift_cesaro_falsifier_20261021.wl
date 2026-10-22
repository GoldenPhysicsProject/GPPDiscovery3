(* GPP Oct21: deterministic fixed-shift Chowla insufficiency falsifier.
   a_n=1 iff n=floor(k^(3/2)); else 0.
   Every FIXED h has only finitely many pairs a_n a_(n+h)=1,
   but the true cumulative Cesaro energy E_N~3N^(1/3).
   Exact proof: gaps s_(k+1)-s_k ->∞,
   A(x)=#{s_k<=x}=x^(2/3)+O(1), integrate A(x)^2/x².
   Numerical checks are diagnostic only. RH OPEN. *)
ClearAll[sp,make,energy,coeff,prefix,nmax,n,h,kmax];
sp[k_Integer]:=Floor[Sqrt[k^3]];
make[m_Integer]:=Module[{kmax,arr},
  kmax=Ceiling[(m+1)^(2/3)]+3;
  arr=ConstantArray[0,m];
  Do[If[1<=sp[k]<=m,arr[[sp[k]]]=1],{k,1,kmax}];
  arr
];
nmax=100000;
coeff=make[nmax];
prefix=Accumulate[coeff];
Print[{"support_count_at_100000",prefix[[-1]],"growth_reference",N[nmax^(2/3),15]}];
Do[
 Module[{val},
  val=N[Sum[prefix[[r]]^2/(r(r+1)),{r,1,n-1}],18];
  Print[{"N",n,"E_N",val,"E_N/N^(1/3)",N[val/n^(1/3),15]}]
 ],{n,{1000,10000,100000}}];
Do[
 Module[{cnt=Total[coeff[[1;;nmax-h]]*coeff[[1+h;;nmax]]]},
   Print[{"shift",h,"pair_count_to_100000",cnt,"fixed-shift-bound-check",cnt<=prefix[[-1]]}]
 ],{h,{1,2,3,5,10,30,100}}];
Print["GENERAL all-N proof uses strict increasing gaps s_(k+1)-s_k >= (k+1)^(3/2)-k^(3/2)-1 ->infinity."];
Print["Every fixed-shift correlation is EVENTUALLY CONSTANT, yet E_N=3N^(1/3)+O(1)."];
Print["Not a multiplicative sequence and not a counterexample to arithmetic Mobius correlation theorems."];
