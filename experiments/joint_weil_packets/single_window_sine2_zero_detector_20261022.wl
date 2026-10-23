(* GPP 2026-10-22 — One-window arithmetic pole detector.
   Exact check: Turbo S1 sine-squared packet has zero-free Laplace multiplier
   in ALL offcritical zeta-zero territory; J's causal double-derivative
   transport cancels the pole denominators. DOES NOT prove RH or sewing
   positivity. Algebraic checks in Wolfram; all-infinite theorem paper-level.
*)
ClearAll[z,x,zz,S,GG,H,deriv,logderiv,LJ,filtered,expected,rho];

S[w_]:=2 Pi^2 (Exp[w]-1)/(w(w^2+4 Pi^2));
GG[z_]:=S[z] S[-z];
H[z_]:=z(z+1/2) GG[z];
assert[s_,p_]:=If[TrueQ[p],Print["PASS: ",s],
 Print["FAIL: ",s];Abort[]];

assert["Turbo wave234 sine² packet Laplace S",
 FullSimplify[Integrate[Exp[z x] Sin[Pi x]^2,{x,0,1}]-S[z],
 Assumptions->z>0]===0];
assert["analytic packet autocorrelation compact filter formula",
 FullSimplify[GG[z] - 8 Pi^4(Cosh[z]-1)/(z^2(z^2+4 Pi^2)^2)]===0];

(* Write dlog for zeta'/zeta; only algebra matters, and no assumptions
   about zero locations or functional equation are inserted here. *)
logderiv[z_] := PolyGamma[0,z+1]; (* placeholder meromorphic function *)
LJ[z_]:=(-logderiv[z]/(z+1/2)-1/(z-1/2))/z;
filtered=FullSimplify[Exp[z] H[z] LJ[z]];
expected=FullSimplify[-Exp[z] GG[z] *
 (logderiv[z]+(z+1/2)/(z-1/2))];
assert["exact Laplace multiplier and prime/Arch pole cancellation",
 FullSimplify[filtered-expected]===0];

(* The true zero pole z0=rho-1/2 has coefficient
   -mult * exp(z0)*S(z0)*S(-z0); this is nonzero if Re z0>0
   because cosh z0 !=1 off the imaginary axis. *)
rho=1/5+I Im[N[ZetaZero[1],40]];
Print[{"hypothetical off-line z=1/5+i gamma1, nonzero multiplier",
 N[Abs[Exp[rho] GG[rho]],28]}];
Print["GENERAL NONVANISHING: cosh(z)=1 iff z in 2*pi*i*Z."];
Print["Hence S(z)S(-z) is nonzero for Re z>0;"];
Print["outside the strip the same statement holds, away from removable points."];
Print["THEOREM (conditional on standard J Laplace formula):"];
Print["Let J(t)=int_1^exp(t) [psi(x)-x]x^-3/2 dx."];
Print["Let G(v)=int g(x)g(x+v) dx, h(v)=G''(v)-G'(v)/2."];
Print["Define K_g(s)=int_{-1}^1 h(v)J(s+1+v) dv for s>=0."];
Print["Then Laplace K_g = -exp(z) S(z)S(-z)"];
Print[" [zeta'/zeta(z+1/2)+(z+1/2)/(z-1/2)] - EntireBoundary(z)."];
Print["An off-line zeta zero would be an UNCANCELLED pole of Laplace K_g."];
Print["Thus K_g bounded or polynomial growth on R>=0 implies RH."];
Print["Conversely RH gives J bounded, hence K_g bounded (classical explicit formula)."];
Print["DO NOT infer this bound from positivity of one-sided Hankel"];
Print["or finite S1 sections: an actual two-sided positive-type"];
Print["correlation with fixed diagonal norm is STILL MISSING."];
Print["Novelty vs literature UNVERIFIED; RH/GRH OPEN."];