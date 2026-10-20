"""Uniform finite-block moment-tail replacement bound; floating evaluation.
No interval certificate and no control of Fourier modes outside this block.
"""
import json
import mpmath as mp
from packet_cross import cross
from modulated_interval import events
mp.mp.dps=100
L=mp.mpf(4);ns=list(range(-8,9));a=[j*mp.pi/L for j in ns];A=max(abs(t) for t in a);ev=events(L)
s=mp.fsum(w/(1-u/(2*L)) for u,w in ev)
b=abs(mp.digamma(mp.mpf('.25'))-mp.log(mp.pi))+4+mp.log(3)+s
Q=mp.matrix([[cross(L,t,u,ev) for u in a] for t in a]);lam=mp.eigsy(Q,eigvals_only=True)[0]
rows=[]
for Ts in ['100','1000','5000']:
 T=mp.mpf(Ts)
 def I(p):return T**(1-p)*((mp.log(T)+b)/(p-1)+1/mp.mpf(p-1)**2)
 for m in [8,16,24,32,48]:
  rho=4/L*mp.fsum(abs(t)**(2*m) for t in a)*I(2*m+2)/(1-A/T)**2
  # trace(D* G D), with G_rq=4/L I_(r+q+2) for even r+q, zero otherwise
  trace=4/L*mp.fsum(mp.fsum(t**(r+q) for t in a)*I(r+q+2) for r in range(m) for q in range(m) if (r+q)%2==0)
  error=(2*mp.sqrt(rho*trace)+rho)/(2*mp.pi)
  rows.append({'T':Ts,'m':m,'rho':mp.nstr(rho,18),'trace_S':mp.nstr(trace,18),'uniform_operator_error':mp.nstr(error,18),'error_over_finite_min':mp.nstr(error/lam,12)})
print(json.dumps({'scope':'Uniform on 17-mode block only; analytic inequality, non-certified constants; signed moment-tail matrix not assembled','L':'4','finite_min':mp.nstr(lam,24),'rows':rows},indent=2))
