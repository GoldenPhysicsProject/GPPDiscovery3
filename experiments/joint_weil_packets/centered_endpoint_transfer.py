"""Independent finite arithmetic quadrature of centered endpoint energy identity.
No infinite arithmetic bound, no interval certificate, no Lean certification.
"""
import numpy as np, math, json
from fractions import Fraction as F
# Exact polynomial identity: r=e^(-t/2)P, x=e^(-t/2)Q, Q'=P.
def mul(a,b):
 out=[F(0)]*(len(a)+len(b)-1)
 for i,u in enumerate(a):
  for j,v in enumerate(b):out[i+j]+=u*v
 return out
def sub(a,b):
 return [(a[i] if i<len(a) else F(0))-(b[i] if i<len(b) else F(0)) for i in range(max(len(a),len(b)))]
exact=0
for P in [[F(1)],[F(1),F(-2)],[F(0),F(2),F(3)],[F(1,3),F(-4,7),F(2,5),F(9,2)]]:
 Q=[F(0)]+[v/F(i+1) for i,v in enumerate(P)]
 q2=mul(Q,Q);der=[(i+1)*q2[i+1] for i in range(len(q2)-1)]
 lhs=sub(mul(P,P),mul(sub(P,Q),sub(P,Q)))
 rhs=sub(der,q2)
 assert all(x==0 for x in sub(lhs,rhs));exact+=1
N=100000;lam=np.zeros(N+1);sieve=np.ones(N+1,dtype=bool);sieve[:2]=False
for p in range(2,N+1):
 if sieve[p]:
  sieve[p*p::p]=False;n=p
  while n<=N:lam[n]=math.log(p);n*=p
nn=np.arange(1,N+1,dtype=float)
psi=np.cumsum(lam)[1:];S=np.cumsum(lam[1:]/np.sqrt(nn));A=np.cumsum(lam[1:]/nn)
c=-2*(1+0.577215664901532860606512090082402431)
g,w=np.polynomial.legendre.leggauss(8);rows=[]
for X in [10,100,1000,10000,100000]:
 n=nn[:X-1];lo=np.log(n);hi=np.log(n+1);mid=(lo+hi)/2;half=(hi-lo)/2
 t=mid[:,None]+half[:,None]*g[None,:];h=np.exp(-t/2);root=1/h
 j0=2*S[:X-1,None]-2*h*psi[:X-1,None]-2*root+2
 jhalf=2*A[:X-1,None]-2*h*S[:X-1,None]-2*t+2-2*h
 r=root*(jhalf-c);v=j0-c*h
 Er=float(np.sum(half[:,None]*w[None,:]*r*r));Ev=float(np.sum(half[:,None]*w[None,:]*v*v))
 T=math.log(X);hT=X**(-.5)
 j0T=2*S[X-1]-2*hT*psi[X-1]-2/hT+2
 jhalfT=2*A[X-1]-2*hT*S[X-1]-2*T+2-2*hT
 rT=(jhalfT-c)/hT;x_formula=rT-j0T+c*hT
 xT=float(np.sum(half[:,None]*w[None,:]*np.exp(-(T-t)/2)*r))
 assert abs(xT-x_formula)<1e-8
 error=Er-Ev-xT*xT
 assert abs(error)<1e-7
 rows.append({'X':X,'weighted_endpoint_energy':Er,'centered_J0_energy':Ev,'boundary_square':xT*xT,'identity_error':error,'r_at_T':rT})
print(json.dumps({'scope':__doc__,'exact_polynomial_checks':exact,'rows':rows},indent=2))
