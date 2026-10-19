"""Parity-resolved shifted-product control. Numerical, not certified."""
import json,mpmath as mp
from packet_cross import cross,arch_cross
from modulated_interval import events
def fake(L,t,s,ev,theta):
    d=t-s
    if not d: prime=mp.fsum(2*mp.cosh(theta*u)*w*mp.cos(t*u) for u,w in ev)
    else: prime=mp.fsum(2*mp.cosh(theta*u)*w/(1-u/(2*L))*mp.cos((t+s)*u/2)*mp.sin(d*(L-u/2))/(L*d) for u,w in ev)
    val=-prime
    for a in [mp.mpf(1)/2+theta,mp.mpf(1)/2-theta]:
        ell=lambda v:2*mp.sinh((a+1j*v)*L)/(a+1j*v)/mp.sqrt(2*L)
        val+=arch_cross(L,t,s,a)+2*mp.re(ell(t)*ell(s))
    return val

mp.mp.dps=85
rows=[]
for Ls in ['0.5','1.4','2']:
 L=mp.mpf(Ls); n=8;ev=events(L); ts=[j*mp.pi/L for j in range(-n,n+1)]
 E=mp.matrix(17,9);O=mp.matrix(17,8);E[8,0]=1
 for j in range(1,9):
  E[8+j,j]=E[8-j,j]=1/mp.sqrt(2)
  O[8+j,j-1]=1/mp.sqrt(2);O[8-j,j-1]=-1/mp.sqrt(2)
 theta=mp.mpf(1)/4
 Q=mp.matrix([[fake(L,t,s,ev,theta) for s in ts] for t in ts])
 P=mp.matrix(17)
 for a in [mp.mpf('.5')+theta,mp.mpf('.5')-theta]:
  ell=lambda t:2*mp.sinh((a+1j*t)*L)/(a+1j*t)/mp.sqrt(2*L)
  P+=mp.matrix([[2*mp.re(ell(t)*ell(s)) for s in ts] for t in ts])
 row={'L':Ls,'theta':'1/4'}
 for name,V in [('even',E),('odd',O)]:
  R=V.T*(Q-P)*V;W=V.T*Q*V
  er=mp.eigsy(R,eigvals_only=True);eq=mp.eigsy(W,eigvals_only=True)
  row[name]={'negative_R':sum(x<0 for x in er),'negative_Q':sum(x<0 for x in eq),'eigenvalues_R':[mp.nstr(x,18) for x in er],'eigenvalues_Q':[mp.nstr(x,18) for x in eq]}
 rows.append(row)
print(json.dumps({'scope':'85-digit finite 17-mode control, no interval certificate','rows':rows},indent=2))
