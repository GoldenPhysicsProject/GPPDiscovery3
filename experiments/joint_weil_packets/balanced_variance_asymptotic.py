"""Finite double-precision diagnostics, not certificates or RH evidence."""
import numpy as np, math, json
N=10_000_000
sieve=np.ones(N+1,dtype=bool);sieve[:2]=False
for p in range(2,math.isqrt(N)+1):
 if sieve[p]:sieve[p*p::p]=False
primes=np.flatnonzero(sieve);rows=[]
for X in [10000,100000,1000000,10000000]:
 t=math.log(X);pp=primes[primes<=X];lp=np.log(pp);a=2*lp*(1/np.sqrt(pp)-1/math.sqrt(X))
 for eps in [.5,.5-.4/math.sqrt(t)]:
  r=pp.astype(float)**(-eps)
  variance=float(np.sum(a*a*(1-r*r)))
  mean=float(np.sum(a*r))
  for p in pp[pp*pp<=X]:
   p=int(p);r0=p**(-eps);powers=[];n=p
   while n<=X:powers.append(n);n*=p
   coeff=[2*math.log(p)*(n**(-.5)-X**(-.5)) for n in powers]
   mean+=sum(coeff[k]*r0**(k+1) for k in range(1,len(coeff)))
   variance+=sum(coeff[k]*coeff[l]*(r0**abs(k-l)-r0**(k+l+2)) for k in range(len(coeff)) for l in range(len(coeff)) if k+l>0)
  def integ(exponent):return math.expm1(exponent*t)/exponent if exponent else t
  C=2*(integ(.5-eps)-X**(-.5)*integ(1-eps))+2*(1-X**(-.5))
  rows.append({'X':X,'epsilon':eps,'signed_mean':mean-C,'variance':variance,'V_minus_2t2_plus_12t':variance-2*t*t+12*t})
print(json.dumps({'scope':__doc__,'rows':rows},indent=2))
