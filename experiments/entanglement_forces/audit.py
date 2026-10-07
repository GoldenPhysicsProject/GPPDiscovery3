"""Exact finite arithmetic and floating-point information-geometry checks.

No physical gauge or gravity identification is assumed; standard library only.
"""
from fractions import Fraction as F
from collections import Counter
from math import log, sin, cos
import json

def dimensions(weights):
    mult = list(Counter(weights).values())
    return {'multiplicities': mult, 'unitary_commutant_dim': sum(m*m for m in mult),
            'special_unitary_commutant_dim': sum(m*m for m in mult)-1}

spectra = {
    'prime_2_first_four_unnormalized': [F(1,2)**n for n in range(4)],
    'three_plus_one': [F(1,5)]*3+[F(2,5)],
    'two_plus_two': [F(1,6)]*2+[F(1,3)]*2,
    'four_equal': [F(1,4)]*4,
}
out = {'stabilizers': {k:dimensions(v) for k,v in spectra.items()}}
assert [out['stabilizers'][k]['special_unitary_commutant_dim'] for k in spectra] == [3,9,7,15]

# Exact product-energy shells log(m)+log(n)=log(k).
shells = {}
for k in (4,6,9,25,30):
    states=[(m,k//m) for m in range(1,k+1) if k%m==0]
    assert all(m*n==k for m,n in states)
    ratios=[F(m,n) for m,n in states]
    assert len(set(ratios))==len(states)  # difference energies split the shell
    shells[str(k)]={'states':states,'dimension':len(states)}
assert shells['25']['dimension']==3
out['joint_energy_shells']=shells

# Exact relative entropy for rotation in a two-dimensional eigenspace.
rot=[]
for a,b in ((.2,.2),(.2,.4),(1/3,1/7)):
    theta=.37
    c2,s2=cos(theta)**2,sin(theta)**2
    direct=a*log(a)+b*log(b)-(c2*a+s2*b)*log(a)-(s2*a+c2*b)*log(b)
    formula=s2*(a-b)*log(a/b)
    assert abs(direct-formula)<1e-14
    assert formula>=0
    rot.append({'a':a,'b':b,'relative_entropy':formula})
out['rotation_cost']=rot

# First-order response of two initially independent finite prime occupations
# under probabilities proportional to exp(-epsilon * n*m).
def moments(p,q,eps):
    from math import exp
    rows=[(n,m,p**(-n)*q**(-m)*exp(-eps*n*m)) for n in range(12) for m in range(12)]
    Z=sum(w for n,m,w in rows)
    en=sum(n*w for n,m,w in rows)/Z
    em=sum(m*w for n,m,w in rows)/Z
    cov=sum(n*m*w for n,m,w in rows)/Z-en*em
    vn=sum(n*n*w for n,m,w in rows)/Z-en*en
    vm=sum(m*m*w for n,m,w in rows)/Z-em*em
    return cov,vn,vm
h=1e-5
c,vn,vm=moments(2,3,0)
deriv=(moments(2,3,h)[0]-moments(2,3,-h)[0])/(2*h)
assert abs(c)<1e-14
assert abs(deriv+vn*vm)<1e-7
out['interaction_response']={'finite_cutoff':11,'cov_at_zero':c,'derivative':deriv,'exact_derivative_for_truncation':-vn*vm}
print(json.dumps(out,indent=2))
