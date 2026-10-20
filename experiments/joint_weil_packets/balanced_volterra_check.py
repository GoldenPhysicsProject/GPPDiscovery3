"""Independent finite arithmetic checks of balanced regulator identity. No RH bound."""
import math,json
import mpmath as mp
mp.mp.dps=60
N=1000;lam={}
for p in range(2,N+1):
 if all(p%d for d in range(2,math.isqrt(p)+1)):
  n=p
  while n<=N:lam[n]=mp.log(p);n*=p
rows=[]
def phi(a,d):return d if not a else mp.expm1(a*d)/a
for X in [10,100,1000]:
 T=mp.log(X)
 for eps in map(mp.mpf,['0','.01','.25','.5','1']):
  M=2*mp.fsum(w*mp.power(n,-eps)*(1/mp.sqrt(n)-1/mp.sqrt(X)) for n,w in lam.items() if n<=X)
  C=2*(phi(mp.mpf('.5')-eps,T)-mp.exp(-T/2)*phi(1-eps,T))+2*(1-mp.exp(-T/2))
  direct=M-C
  points=sorted({1,X}|{n for n in lam if n<X});P=mp.mpf(0);W=mp.mpf(0);i0=mp.mpf(0);i1=mp.mpf(0)
  for a,b in zip(points,points[1:]):
   w=lam.get(a,0);P+=w;W+=w/mp.sqrt(a);A=2*W+2;B=-2*P
   lo=mp.log(a);d=mp.log(mp.mpf(b)/a)
   def cell(q):return mp.exp(q*lo)*phi(q,d)
   i0+=A*cell(-eps)+B*cell(-eps-mp.mpf('.5'))-2*cell(mp.mpf('.5')-eps)
   i1+=A*cell(mp.mpf('.5')-eps)+B*cell(-eps)-2*cell(1-eps)
  J=2*mp.fsum(w*(1/mp.sqrt(n)-1/mp.sqrt(X)) for n,w in lam.items() if n<=X)-2*mp.sqrt(X)+2
  via=mp.exp(-eps*T)*J+eps*(1+2*eps)*i0+eps*(1-2*eps)*mp.exp(-T/2)*i1
  err=abs(direct-via);assert err<mp.mpf('1e-50')
  rows.append(dict(X=X,epsilon=str(eps),balanced=mp.nstr(direct,25),volterra_error=mp.nstr(err,6)))
print(json.dumps(dict(scope='Finite analytic-cell identity checks, not uniform estimates or interval/Lean certification',rows=rows),indent=2))
