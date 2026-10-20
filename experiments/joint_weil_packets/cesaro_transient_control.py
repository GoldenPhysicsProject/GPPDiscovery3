"""Remove explicit deterministic transient; compare finite zero list independently."""
import json,math,struct
from pathlib import Path
import mpmath as mp
mp.mp.dps=60
N=100000;sieve=bytearray(b'\1')*(N+1);sieve[:2]=b'\0\0'
for p in range(2,math.isqrt(N)+1):
 if sieve[p]:sieve[p*p:N+1:p]=b'\0'*((N-p*p)//p+1)
lam={}
for p in range(2,N+1):
 if sieve[p]:
  n=p
  while n<=N:lam[n]=mp.log(p);n*=p
zeros=[]
for name in ['zeros_14.dat','zeros_5000.dat']:
 raw=Path(__file__).with_name(name).read_bytes();pos=8
 for k in range(struct.unpack_from('<Q',raw)[0]):
  t0,t1,n0,n1=struct.unpack_from('<ddQQ',raw,pos);pos+=32;z=0
  for j in range(n1-n0):
   z+=int.from_bytes(raw[pos:pos+13],'little');pos+=13
   zeros.append(mp.mpf(t0)+mp.mpf(z)*mp.mpf(2)**-101)
 assert pos==len(raw)
kappa=2+mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi)
def transient(t):
 q=mp.exp(-2*t)
 # Exact identity sum q^n/[n(4n+1)] = -log(1-q) - 4*q/5 * 2F1(1,5/4;9/4;q).
 S=-mp.log1p(-q)-4*q/5*mp.hyp2f1(1,mp.mpf(5)/4,mp.mpf(9)/4,q)
 return mp.exp(-t/2)*(2*mp.log(2*mp.pi)-S)
rows=[]
for X in [10,100,1000,10000,100000]:
 t=mp.log(X);J=2*mp.fsum(w*(1/mp.sqrt(n)-1/mp.sqrt(X)) for n,w in lam.items() if n<=X)-2*mp.sqrt(X)+2
 R=transient(t);B=J-kappa-R
 derivative=mp.diff(transient,t)
 expected=-mp.exp(-t/2)*(mp.log(2*mp.pi)+mp.log1p(-mp.exp(-2*t))/2)
 assert abs(derivative-expected)<mp.mpf('1e-50')
 terms=[-2*mp.re(mp.exp(1j*g*t)/((mp.mpf('.5')+1j*g)*(1j*g))) for g in zeros]
 small=mp.fsum(terms[:4520]);large=mp.fsum(terms)
 rows.append(dict(X=X,J=mp.nstr(J,25),transient=mp.nstr(R,25),arithmetic_residual=mp.nstr(B,25),listed_4520=mp.nstr(small,25),listed_30324=mp.nstr(large,25),error_4520=mp.nstr(abs(B-small),10),error_30324=mp.nstr(abs(B-large),10)))
print(json.dumps(dict(scope='Prime-defined residual compared to finite listed on-line zeros; no omitted-zero bound, RH or interval certificate',kappa=mp.nstr(kappa,30),zero_count=len(zeros),rows=rows),indent=2))
