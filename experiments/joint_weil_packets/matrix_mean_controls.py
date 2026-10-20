"""Noncommuting mean controls. Exact rational first example; floating radical second."""
from fractions import Fraction as F
import json
import mpmath as mp
mp.mp.dps=70
# Paper Def 3.2, n=1: P0=diag(2,1),P1=I,S0=ones,S1=0.
left=[[F(1,3),F(1,3)],[F(1,2),F(1,2)]]
assert left[0][1]!=left[1][0]
# Sandwich repair with alternating P=diag(2,1),Q=diag(1,2).
# For a CONSTANT X=ones and any even number of weights, output Y below.
x=mp.matrix([[1,1],[1,1]]);b=2*mp.sqrt(2)/3
y=mp.matrix([[1,b],[b,1]])
ev=mp.eigsy(y,eigvals_only=True)
assert min(ev)>0
assert abs(b-1)>mp.mpf('.05')
# Scalar Cesaro plus an explicit error envelope: Xk=X+(-1)^k/(k+1) I.
rows=[]
for N in [10,100,1000]:
 avgerr=mp.fsum((-1)**k/mp.mpf(k+1) for k in range(N))/N
 envelope=mp.fsum(1/mp.mpf(k+1) for k in range(N))/N
 assert abs(avgerr)<=envelope
 rows.append(dict(N=N,actual_operator_error=mp.nstr(abs(avgerr),20),propagated_error_bound=mp.nstr(envelope,20)))
print(json.dumps(dict(scope='Model tests; no RH or actual arithmetic positivity claim',left_weighted_mean=[[str(z) for z in row] for row in left],sandwich_constant_offdiagonal=mp.nstr(b,30),desired_constant_offdiagonal='1',sandwich_eigenvalues=[mp.nstr(z,30) for z in ev],scalar_controls=rows),indent=2))
