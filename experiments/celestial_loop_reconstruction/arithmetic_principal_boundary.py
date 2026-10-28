"""Positive Archimedean-weighted principal-series CAR carrier; no interacting gravity."""
import cmath,math,json
import numpy as np
from scipy.special import digamma
mu=1.
def kinetic(lam):return float((digamma(.25+.5j*lam)-digamma(.25)).real)
def weight(lam):return 1/(mu+kinetic(lam))
def transport(g,z,lam,k):
 a,b,c,d=g;D=c*z+d
 return (a*z+b)/D,abs(D)**(2+2j*lam)*(D/abs(D))**k,abs(D)**(-4)
gs=[(1,.4+.2j,0,1),(2,0,0,.5),(1,0,.3+.1j,1),(0,-1,1,0)]
res=0
for g in gs:
 for z in [.2+.4j,1.2-.7j,-.8+.9j]:
  for lam in [-3.,0.,2.7]:
   for k in [0,1,-1]:
    zp,A,jac=transport(g,z,lam,k)
    res=max(res,abs(abs(A)**2*jac-1))
    assert abs(abs(A)**2*jac-1)<1e-12
    _,An,_=transport(tuple(-x for x in g),z,lam,k)
    assert abs(An-(-1)**k*A)<1e-12
# Coupled composition transports: g then h. Spin lift retains the central sign.
def product(h,g):
 a,b,c,d=h;e,f,j,k=g
 return (a*e+b*j,a*f+b*k,c*e+d*j,c*f+d*k)
for g in gs:
 for h in gs:
  z=.37+.61j;lam=.8
  zg,A,_=transport(g,z,lam,1);_,B,_=transport(h,zg,lam,1)
  _,C,_=transport(product(h,g),z,lam,1)
  assert abs(A*B-C)<1e-11
for lam in np.linspace(-20,20,401):
 assert kinetic(lam)>-1e-12 and 0<weight(lam)<=1+1e-12
 assert abs(weight(lam)-weight(-lam))<1e-13
# Weighted finite packet Gram and fermionic Wick four-point determinant.
rng=np.random.default_rng(17);F=rng.normal(size=(41,4))+1j*rng.normal(size=(41,4))
w=np.array([weight(x) for x in np.linspace(-5,5,41)])
G=F.conj().T@(w[:,None]*F)
assert np.linalg.eigvalsh(G).min()>0
four=G[0,2]*G[1,3]-G[0,3]*G[1,2]
connected=four-(G[0,2]*G[1,3]-G[0,3]*G[1,2])
assert connected==0
result=dict(mobius_density_residual=res,composition='passed',spin_half_central_sign='passed',positive_even_archimedean_weight='passed',minimum_packet_gram_eigenvalue=float(np.linalg.eigvalsh(G).min()),connected_four_point=[connected.real,connected.imag],scope='Principal-series Gaussian CAR boundary carrier only; primes not coupled; no translations, local CFT, or interacting gravitational reconstruction')
print(json.dumps(result,indent=2));open('arithmetic_principal_boundary_results.json','w').write(json.dumps(result,indent=2)+'\n')
