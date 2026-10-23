"""Exact finite-rank KMS checks and local OS Gram checks. Not full Weil positivity."""
from fractions import Fraction as F
import mpmath as m
import json
m.mp.dps=50
count=0
for p in [2,3,5]:
 d=5
 for i in range(d):
  for j in range(d):
   # A=E_ij, B=E_ji: omega(A alpha_i(B))=omega(BA).
   rho_i=(1-F(1,p))*F(1,p)**i
   rho_j=(1-F(1,p))*F(1,p)**j
   assert rho_i*F(p)**(i-j)==rho_j
   count+=1
 # Local OS: arbitrary real finite-rank observables at times in [0, beta/2].
 times=[m.mpf('0'),m.mpf('.125'),m.mpf('.25'),m.mpf('.5')]
 As=[m.matrix(d,d) for _ in times]
 for a,A in enumerate(As):
  for i in range(d):
   for j in range(d): A[i,j]=((i+2*j+a)%5)-2
 X=[]
 for tau,A in zip(times,As):
  X.append(m.matrix([[m.power(p,-tau*i)*A[i,j]*m.power(p,-(m.mpf('.5')-tau)*j)*m.sqrt(1-m.mpf(1)/p) for j in range(d)] for i in range(d)]))
 G=m.matrix([[sum(X[a][i,j]*X[b][i,j] for i in range(d) for j in range(d)) for b in range(4)] for a in range(4)])
 ev=m.eigsy(G,eigvals_only=True)
 assert ev[0]>0
 print(json.dumps(dict(p=p,min_OS_eigenvalue=str(ev[0]))))
print(json.dumps(dict(exact_KMS_matrix_unit_checks=count,scope=__doc__)))
