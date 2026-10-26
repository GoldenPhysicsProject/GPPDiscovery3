"""Exact rational audit of multiple-zero cardinal interpolation. Not an RH test."""
from fractions import Fraction as Q
nodes = (Q(1,4), Q(-1,4), Q(0))
multiplicity = (2,2,1)
def cardinal(j,x):
    result=Q(1)
    for k,z in enumerate(nodes):
        if k != j:
            result *= ((x-z)/(nodes[j]-z))**multiplicity[k]
    return result
values=[cardinal(0,x)-cardinal(1,x) for x in nodes]
assert values == [1,-1,0]
form=2*values[0]*values[1]+2*values[1]*values[0]+values[2]**2
assert form == -4
# Positivity on a transform subspace does not imply truncated positivity.
# Every (a,-a,2a) has full form 2*a*a >=0, but pair form -2*a*a.
for a in map(Q,range(-10,11)):
    assert 2*a*(-a)+(2*a)**2 >= 0
    if a:
        assert 2*a*(-a)<0
print("PASS: multiplicity-two cardinal values",values,"Q =",form)
print("PASS: positive full restricted form can have negative pair truncation")
