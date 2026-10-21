"""Exact rational tests of finite-regulator prime-power covariance factorization."""
from fractions import Fraction as F
import json
rows=[]
for r in [F(0),F(1,10),F(1,2),F(9,10),F(1)]:
 for a in [[F(2),F(3)], [F(1),F(-2),F(3),F(-4)], [F(1,3),F(2,7),F(5,11)]]:
  n=len(a)
  cov=sum(a[k-1]*a[l-1]*(r**abs(k-l)-r**(k+l)) for k in range(1,n+1) for l in range(1,n+1))
  gram=(1-r*r)*sum(sum(a[k-1]*r**(k-j) for k in range(j,n+1))**2 for j in range(1,n+1))
  assert cov==gram and gram>=0
  if all(x>=0 for x in a):
   slope=sum(sum(a[k-1] for k in range(j,n+1))**2 for j in range(1,n+1))
   assert cov<=(1-r*r)*slope
  rows.append({'r':str(r),'coefficients':list(map(str,a)),'variance':str(cov),'identity_exact':True})
print(json.dumps({'scope':'Exact rational finite identities; no arithmetic mean bound', 'cases':rows},indent=2))
