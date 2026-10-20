"""Finite prime identity and step-integral check, no asymptotic/RH claim."""
import math,json
import mpmath as mp
mp.mp.dps=60
rows=[]
for Xs in ['2','10','100','1000']:
 X=mp.mpf(Xs);N=int(X); lam={}
 for p in range(2,N+1):
  if all(p%d for d in range(2,math.isqrt(p)+1)):
   n=p
   while n<=N:lam[n]=mp.log(p);n*=p
 psi=mp.fsum(lam.values())
 J=2*mp.fsum(w/mp.sqrt(n) for n,w in lam.items())-2*psi/mp.sqrt(X)-2*mp.sqrt(X)+2
 points=sorted(set([mp.mpf(1),X]+[mp.mpf(n) for n in lam]))
 direct=mp.mpf(0);current=mp.mpf(0)
 for a,b in zip(points,points[1:]):
  current+=lam.get(int(a),0)
  direct+=2*current*(1/mp.sqrt(a)-1/mp.sqrt(b))-2*(mp.sqrt(b)-mp.sqrt(a))
 assert abs(J-direct)<mp.mpf('1e-55')
 # exact finite divisor identity Lambda(n)=sum_{d|n}mu(d)log(n/d)
 def mu(n):
  v=1;d=2
  while d*d<=n:
   if n%d==0:
    n//=d;v=-v
    if n%d==0:return 0
   d+=1
  return -v if n>1 else v
 conv=2*mp.fsum(mu(d)*mp.fsum(mp.log(m)*(1/mp.sqrt(d*m)-1/mp.sqrt(X)) for m in range(1,N//d+1)) for d in range(1,N+1))-2*mp.sqrt(X)+2
 assert abs(J-conv)<mp.mpf('1e-55')
 rows.append({'X':Xs,'J':mp.nstr(J,30),'S':mp.nstr(J/mp.log(X),30),'step_error':mp.nstr(abs(J-direct),5),'mobius_error':mp.nstr(abs(J-conv),5)})
print(json.dumps({'scope':'Finite 60-digit numerical identity checks, not interval certified','rows':rows},indent=2))
