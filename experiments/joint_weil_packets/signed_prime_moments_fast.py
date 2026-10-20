"""Signed prime moments. Recurrence with guard digits; NOT interval certified."""
import json,time
import mpmath as mp
from modulated_interval import events
mp.mp.dps=260
L=mp.mpf(4); T=mp.mpf(1000)
primes=[(u,w/(1-u/(2*L))) for u,w in events(L)]
def moments(x):
 x=abs(x)
 if not x:return {q:1/mp.mpf(q-1) for q in range(2,49,2)}
 e=mp.exp(1j*x)
 E=-mp.ci(x)+1j*(mp.pi/2-mp.si(x))
 out={}
 for q in range(1,48):
  E=(e+1j*x*E)/q
  if (q+1)%2==0:out[q+1]=mp.re(E)
 return out
start=time.time(); sums={q:[] for q in range(2,49,2)}
for u,c in primes:
 a,b,d=[moments(x*T) for x in (u,u+2*L,u-2*L)]
 for q in sums:sums[q].append(c*(a[q]/2-b[q]/4-d[q]/4))
checks=[]
for x in [mp.mpf('.3'),mp.mpf('1000'),mp.mpf('15999')]:
 vals=moments(x)
 for q in [2,24,48]:
  error=abs(vals[q]-mp.re(mp.expint(q,-1j*x)))
  assert error<mp.mpf('1e-70')
  checks.append(dict(x=str(x),q=q,error=mp.nstr(error,8)))
rows=[]
for q,terms in sums.items():
 v=mp.fsum(terms); env=mp.fsum(c for u,c in primes)/(q-1)
 rows.append(dict(q=q,normalized_signed_prime_moment=mp.nstr(v,50),prime_moment=mp.nstr(v*T**(1-q),50),absolute_envelope=mp.nstr(env*T**(1-q),20),signed_to_envelope=mp.nstr(v/env,15)))
print(json.dumps(dict(scope='Floating signed PRIME contribution only, to be SUBTRACTED from Arch; no PSD or RH certificate',L='4',T='1000',prime_power_count=len(primes),dps=mp.mp.dps,seconds=time.time()-start,checks=checks,rows=rows),indent=2))
