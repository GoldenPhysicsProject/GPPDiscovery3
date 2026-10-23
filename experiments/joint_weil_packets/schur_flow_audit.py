"""Exact rational tests: growing-window inverse and directional energy."""
from fractions import Fraction as F
import json
checks=0
rows=[]
for c in [F(1,3),F(1,2),F(2,3)]:
 for L in [F(1,5),F(1,2),F(3,4)]:
  for h in [F(1,100),F(1,20)]:
   assert 1-c*(L+h)>0
   old=c/(1-c*L);new=c/(1-c*(L+h))
   exact=c*c*h/((1-c*L)*(1-c*(L+h)))
   assert new-old==exact
   s=1-c*h/(1-c*L)
   # Integral operator I-c|1><1|; normalized new-cell coupling.
   assert s==(1-c*(L+h))/(1-c*L)
   assert exact==c*c*h/((1-c*L)**2*s)
   checks+=3
   rows.append(dict(c=str(c),L=str(L),h=str(h),inverse_increment=str(exact),schur=str(s)))
# Noncoercive diagonal A_jj=1/j^2. b_j=1/j^2.
directional=[]
for n in [10,100,1000]:
 q=sum((F(1,j*j) for j in range(1,n+1)),F(0))
 # elementary tail sum_{j>n}1/j^2 <= integral_n^infinity x^-2 dx
 assert q<2
 directional.append(dict(n=n,q=float(q),tail_upper_bound=float(F(1,n)),inverse_vector_norm_squared=n))
print(json.dumps(dict(scope=__doc__,exact_rational_checks=checks,rank_one_cases=rows,directional_energy=directional),indent=2))
