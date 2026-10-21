"""Prime-side diagnostic of the coherent mean; no zero input or RH assumption."""
import math,json
import mpmath as mp
mp.mp.dps=50
N=100000
s=bytearray(b'\1')*(N+1);s[:2]=b'\0\0'
for p in range(2,math.isqrt(N)+1):
 if s[p]:s[p*p:N+1:p]=b'\0'*((N-p*p)//p+1)
primes=[p for p in range(2,N+1) if s[p]]
rows=[]
for X in [10,100,1000,10000,100000]:
 t=mp.log(X); root=mp.sqrt(X); layers=[mp.mpf(0)]*3
 for p in primes:
  if p>X:break
  w=2*mp.log(p);n=p;k=1
  while n<=X:
   layers[min(k-1,2)]+=w*(1/mp.sqrt(n)-1/root)
   n*=p;k+=1
 primary=layers[0]-2*root+2
 J=primary+layers[1]+layers[2]
 assert 0<=layers[2]<40
 rows.append(dict(X=X,t=mp.nstr(t,20),primary_signed=mp.nstr(primary,25),
                  prime_square=mp.nstr(layers[1],25),higher_powers=mp.nstr(layers[2],25),
                  coherent_mean=mp.nstr(J,25),square_over_t=mp.nstr(layers[1]/t,15)))
print(json.dumps(dict(scope='50-digit finite arithmetic diagnostics, not certified bounds or a proof of RH',rows=rows),indent=2))
