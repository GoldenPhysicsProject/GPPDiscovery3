"""Endpoint-moment tail enclosure formulas, floating evaluation only.
Weight envelope: |w(t)| <= log(|t|)+b for |t|>=T>=1.
Outputs integrate over both tails WITHOUT the Plancherel factor 1/(2*pi).
"""
import json
import mpmath as mp
from packet_cross import cross
from modulated_interval import events
mp.mp.dps=85
L=mp.mpf(4);ev=events(L);ns=list(range(-8,9));a=[j*mp.pi/L for j in ns]
Q=mp.matrix([[cross(L,t,u,ev) for u in a] for t in a]);vals,V=mp.eigsy(Q);v=V[:,0]
s=mp.fsum(w/(1-u/(2*L)) for u,w in ev)
b=abs(mp.digamma(mp.mpf('.25'))-mp.log(mp.pi))+4+mp.log(3)+s
amax=max(abs(x) for x in a)
M=[mp.fsum((-1)**j*v[k]*a[k]**r for k,j in enumerate(ns)) for r in range(16)]
rows=[]
for Ts in ['100','1000','5000']:
 T=mp.mpf(Ts)
 def I(p):return T**(1-p)*((mp.log(T)+b)/(p-1)+1/mp.mpf(p-1)**2)
 crude=4/L*mp.fsum(abs(x) for x in v)**2/(1-amax/T)**2*I(2)
 for m in [1,4,8,16]:
  B=mp.fsum(abs(v[k])*abs(a[k])**m for k in range(17))
  ER=4/L*B**2/(1-amax/T)**2*I(2*m+2)
  EM=4/L*mp.fsum(abs(M[r])*abs(M[q])*I(r+q+2) for r in range(m) for q in range(m))
  error=2*mp.sqrt(EM*ER)+ER
  rows.append({'T':Ts,'m':m,'moment_tail_upper':mp.nstr(EM,18),'remainder_energy_upper':mp.nstr(ER,18),'quadratic_error_upper':mp.nstr(error,18),'full_tail_upper':mp.nstr(EM+error,18),'crude_tail_upper':mp.nstr(crude,18),'error_over_Q':mp.nstr(error/(2*mp.pi*vals[0]),12)})
print(json.dumps({'scope':'Analytic inequalities with floating constants/eigenvector; no interval certificate, no full-window positivity','L':'4','prime_Q':mp.nstr(vals[0],24),'envelope_b':mp.nstr(b,24),'moments':[mp.nstr(z,18) for z in M],'rows':rows},indent=2))
