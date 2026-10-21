"""Exact coefficient audit of the Selberg convolution identity.
Log primes are treated as independent formal variables; no floating proof.
"""
from collections import defaultdict
import json
N=500
spf=list(range(N+1))
for p in range(2,N+1):
 if spf[p]==p:
  for n in range(p*p,N+1,p):
   if spf[n]==n:spf[n]=p

def fac(n):
 out={}
 while n>1:
  p=spf[n];out[p]=out.get(p,0)+1;n//=p
 return out
factors=[{}]+[fac(n) for n in range(1,N+1)]
def mu(n):
 f=factors[n]
 return 0 if any(v>1 for v in f.values()) else (-1)**len(f)
def square_log(n):
 f=factors[n];out=defaultdict(int)
 for p,a in f.items():
  for q,b in f.items():out[tuple(sorted((p,q)))]+=a*b
 return out
examples=[]
for n in range(1,N+1):
 rhs=defaultdict(int)
 for d in range(1,n+1):
  if n%d==0:
   for key,v in square_log(n//d).items():rhs[key]+=mu(d)*v
 rhs={k:v for k,v in rhs.items() if v}
 f=factors[n];lhs={}
 if len(f)==1:
  p,a=next(iter(f.items()));lhs[(p,p)]=2*a-1
 elif len(f)==2:
  p,q=sorted(f);lhs[(p,q)]=2
 assert rhs==lhs,(n,rhs,lhs)
 assert all(v>=0 for v in rhs.values())
 if n in [8,12,30,49,72]:examples.append({'n':n,'polynomial':{str(k):v for k,v in rhs.items()}})
print(json.dumps({'scope':__doc__,'exact_cases':N,'examples':examples,'infinite_bound_obtained':False},indent=2))
