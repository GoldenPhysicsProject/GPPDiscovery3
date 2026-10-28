from math import isqrt
import json
limit=256
counts=[0]*(limit+1)
for a in range(-isqrt(limit),isqrt(limit)+1):
 for b in range(-isqrt(limit),isqrt(limit)+1):
  n=a*a+b*b
  if 0<n<=limit: counts[n]+=1
chi=lambda n:0 if n%2==0 else (1 if n%4==1 else -1)
for n in range(1,limit+1):
 divisor_sum=sum(chi(d) for d in range(1,n+1) if n%d==0)
 assert counts[n]%4==0
 assert counts[n]//4==divisor_sum
primes=[2,3,5,7,11,13]
rows=[]
for p in primes:
 coeff=[]
 for k in range(5):
  if p**k>limit:break
  actual=counts[p**k]//4
  expected=1 if p==2 else (k+1 if p%4==1 else int(k%2==0))
  assert actual==expected
  coeff.append(actual)
 rows.append({'prime':p,'type':'ramified' if p==2 else 'split' if p%4==1 else 'inert','local_coefficients':coeff})
print(json.dumps({'norms_checked':limit,'identity':'r_2(n)/4=sum_{d|n} chi_-4(d)','prime_checks':rows,'scope':'Exact finite integer checks, not a proof of the scattering formula or analytic continuation.'},indent=2))
