"""Finite independent von Mangoldt and Mobius reconstructions of the signed mean."""
import math,json
import mpmath as mp
mp.mp.dps=50
N=10000
mu=[1]*(N+1);mu[0]=0
prime=[True]*(N+1);prime[0]=prime[1]=False
lam=[mp.mpf(0)]*(N+1)
for p in range(2,N+1):
 if prime[p]:
  for n in range(p,N+1,p):
   mu[n]*=-1
   if n>p:prime[n]=False
  for n in range(p*p,N+1,p*p):mu[n]=0
  n=p
  while n<=N:lam[n]=mp.log(p);n*=p
logs=[mp.mpf(0)]*(N+1); weighted=logs.copy()
for n in range(1,N+1):
 logs[n]=logs[n-1]+mp.log(n)
 weighted[n]=weighted[n-1]+mp.log(n)/mp.sqrt(n)
def kernel(y):
 n=int(mp.floor(y));return 2*(weighted[n]-logs[n]/mp.sqrt(y))
rows=[]
for X in [10,100,1000,10000]:
 x=mp.mpf(X);root=mp.sqrt(x)
 A=mp.fsum(mp.mpf(mu[d])/d*mp.log(x/d) for d in range(1,X+1))
 B=mp.fsum(mp.mpf(mu[d])/d for d in range(1,X+1))
 C=mp.mpf(0);dmin=mp.inf;dmax=-mp.inf
 for d in range(1,X+1):
  y=x/d;D=kernel(y)-2*mp.sqrt(y)*mp.log(y)+6*mp.sqrt(y)
  dmin=min(dmin,D);dmax=max(dmax,D)
  C+=mu[d]*D/mp.sqrt(d)
 direct=2*mp.fsum(lam[n]*(1/mp.sqrt(n)-1/root) for n in range(1,X+1))-2*root+2
 leading=2*root*(A-3*B-1)
 reconstructed=leading+C+2
 assert abs(direct-reconstructed)<mp.mpf('1e-40')
 assert 4<=dmin<=dmax<=12
 rows.append(dict(X=X,leading=mp.nstr(leading,25),mobius_kernel_term=mp.nstr(C,25),
                  J=mp.nstr(direct,25),error=mp.nstr(abs(direct-reconstructed),5),
                  D_min=mp.nstr(dmin,15),D_max=mp.nstr(dmax,15)))
print(json.dumps(dict(scope='50-digit finite identity diagnostics, no RH bound',rows=rows),indent=2))
