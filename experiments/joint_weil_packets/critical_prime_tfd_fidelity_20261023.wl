(* Critical local thermofield product and normalized arithmetic fidelity. *)
ClearAll[e,a,b,z,f];
f[x_,y_]:=Zeta[(x+y)/2]/Sqrt[Zeta[x]Zeta[y]];
Print[FullSimplify[Limit[f[1+e,1+a*e],e->0,Assumptions->a>0],Assumptions->a>0]];
Print[FullSimplify[Limit[f[1+e,b],e->0,Assumptions->b>1],Assumptions->b>1]];
Print[N[f[1+1/10000,1+1/20000],30]];
Print[FullSimplify[(1-z)Sum[z^n,{n,0,Infinity}],Assumptions->0<z<1]];