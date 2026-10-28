"""Exact spinor certificate for four-point interaction ambiguity; not a UV completion."""
from fractions import Fraction as Q
import json,itertools
L=[(Q(1),Q(0)),(Q(0),Q(1)),(Q(1),Q(1)),(Q(1),Q(2))]
R=[(Q(-3),Q(-4)),(Q(-5),Q(-5)),(Q(1),Q(3)),(Q(2),Q(1))]
def br(a,b):return a[0]*b[1]-a[1]*b[0]
def contact(l,r):return br(l[0],l[1])**4*br(r[2],r[3])**4
assert all(sum(L[i][a]*R[i][b] for i in range(4))==0 for a in range(2) for b in range(2))
c=contact(L,R);assert c==625
for i,h in enumerate([-2,-2,2,2]):
 l=L.copy();r=R.copy();t=Q(3,2)
 l[i]=tuple(t*x for x in l[i]);r[i]=tuple(x/t for x in r[i])
 assert contact(l,r)==t**(-2*h)*c
for perm in [(1,0,2,3),(0,1,3,2)]:
 assert contact([L[i] for i in perm],[R[i] for i in perm])==c
# Each leg enters with four spinor powers: uniform soft scaling gives delta^2.
# This polynomial test does not itself impose a momentum-conserving soft path.
for i in range(4):
 l=L.copy();r=R.copy();t=Q(1,3)
 l[i]=tuple(t*x for x in l[i]);r[i]=tuple(t*x for x in r[i])
 assert contact(l,r)==t**4*c
# Multiplying by an arbitrary symmetric polynomial in s,t,u supplies further local terms.
s=br(L[0],L[1])*br(R[1],R[0]);t=br(L[0],L[2])*br(R[2],R[0]);u=br(L[0],L[3])*br(R[3],R[0])
assert s+t+u==0
assert all(sum(x*x for x in p)==s*s+t*t+u*u for p in itertools.permutations([s,t,u]))
result=dict(momentum_conservation='exact',contact=str(c),mandelstam=list(map(str,[s,t,u])),little_group='passed all four legs',identical_helicity_exchanges='passed',soft_polynomial_degree='delta^2 each leg',pole_residues='zero: polynomial has no propagator denominators',scope='Local four-point EFT ambiguity; no UV completion, celestial crossing solution, or loop claim')
print(json.dumps(result,indent=2))
open('four_graviton_contact_audit_results.json','w').write(json.dumps(result,indent=2)+'\n')
