"""Exact rational tests of weighted Schur assembly-error budgets.
Synthetic matrices only; no replay of Grok's dumps and no RH claim.
"""
from fractions import Fraction as Q
import json
# A=diag(10^-80,1), B=(6*10^-40,1)^T, C=B^T A^-1 B+1/2.
a=Q(1,10**80);B=[Q(6,10**40),Q(1)];C=Q(75,2)
base=C-B[0]**2/a-B[1]**2
assert base==Q(1,2)
rows=[]
for soft,hard,e_bound in [(Q(1,10**42),Q(0),Q(1,100)),(-Q(1,10**42),Q(0),Q(1,100)),(Q(1,10**42),Q(1,10000),Q(11,1000)),(Q(1,10**39),Q(0),Q(10))]:
 E=[soft,hard];esq=E[0]**2/a+E[1]**2
 assert esq<=e_bound**2
 new=C-(B[0]+E[0])**2/a-(B[1]+E[1])**2
 loss=2*Q(7)*e_bound+e_bound**2 # b=sqrt(37)<7
 lower=base-loss
 assert new>=lower
 rows.append({'soft_error':str(soft),'hard_error':str(hard),'true_S':str(new),'certified_lower':str(lower),'positive_certificate':lower>0})
assert not rows[-1]['positive_certificate'] and Q(rows[-1]['true_S'])<0
# Conditional S76 budget: b<=7 must be certified separately on the real dump.
nm=200*988;delta=Q(5,10**46);a_cert=Q(2,10**81)
e2=nm*delta**2/a_cert;eupper=Q(1,200)
assert e2<=eupper**2
budget=2*7*eupper+eupper**2
assert budget==Q(2801,40000)
print(json.dumps({'scope':__doc__,'exact_cases':rows,'conditional_s76_example':{'entry_error':str(delta),'e_squared_bound':str(e2),'e_upper':str(eupper),'schur_loss_upper':str(budget),'requires':'certified b<=7, fixed certified A, and separate C error budget'}},indent=2))
