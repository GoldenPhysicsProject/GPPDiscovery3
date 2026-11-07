#!/usr/bin/env python3
"""Exact binary quadratic-form/ideal-class character toy for Q(sqrt(-5)).

Discriminant -20 has class group C2, represented by the two reduced
primitive positive definite binary quadratic forms
  Q0(x,y) = x²+5y², Q1(x,y) = 2x²+2xy+3y².
Since units are ±1, r_Qj(n)/2 counts ideals of norm n in the
corresponding class (a classical arithmetic theorem, NOT proved by
this finite enumeration). Their class-character Fourier channels
a_plus(n)=(r_Q0(n)+r_Q1(n))/2 and
a_minus(n)=(r_Q0(n)-r_Q1(n))/2
correspond to Dedekind zeta and nontrivial class-group Hecke L
respectively. The theorem of global Euler product is mathematical
prior art; the exact coefficient tests here do not establish it.

Under-class-projected Q0's coefficient r_Q0(n)/2 is NOT multiplicative.
It satisfies f(2)=f(3)=0 but f(6)=2. This falsifies naive
"positive integral theta lattice + FE => a single Euler product".
NO global zero calculations or RH proof.
"""
import argparse
from math import gcd, isqrt
import json


def counts(limit: int, cls: int) -> list[int]:
    """Return exact r_Q(n) counts for 1..limit with justified bounds."""
    if limit<6 or cls not in (0,1):
        raise ValueError("limit>=6, cls=0/1 required")
    out=[0]*(limit+1)
    if cls==0:
        xmax=isqrt(limit)
        ymax=isqrt(limit//5)
        for x in range(-xmax,xmax+1):
            for y in range(-ymax,ymax+1):
                n=x*x+5*y*y
                if 1<=n<=limit:
                    out[n]+=1
    else:
        # Q1 = (2x+y)^2/2 + 5y²/2. Thus 5y²<=2N,
        # and |2x+y|<=sqrt(2N); 2x+y has the same parity as y.
        ymax=isqrt(2*limit//5)
        bmax=isqrt(2*limit)
        for y in range(-ymax,ymax+1):
            for u in range(-bmax,bmax+1):
                if (u-y)%2:
                    continue
                x=(u-y)//2
                n=2*x*x+2*x*y+3*y*y
                if 1<=n<=limit:
                    out[n]+=1
    assert all(v%2==0 for v in out)
    return out


def report(limit=128):
    q0=counts(limit,0);q1=counts(limit,1)
    principal=[v//2 for v in q0]
    other=[v//2 for v in q1]
    even=[(q0[n]+q1[n])//2 for n in range(limit+1)]
    odd=[(q0[n]-q1[n])//2 for n in range(limit+1)]
    assert principal[1]==even[1]==odd[1]==1
    assert other[1]==0
    assert (principal[2],principal[3],principal[6])==(0,0,2)
    assert (even[2],even[3],even[6])==(1,2,2)
    assert (odd[2],odd[3],odd[6])==(-1,-2,2)
    assert all(principal[n]+other[n]==even[n] for n in range(1,limit+1))
    assert all(principal[n]-other[n]==odd[n] for n in range(1,limit+1))
    pairs=0
    for a in range(1,limit+1):
        for b in range(1,limit//a+1):
            if gcd(a,b)!=1:
                continue
            for seq in (even,odd):
                assert seq[a*b]==seq[a]*seq[b],(a,b,seq[a*b],seq[a]*seq[b])
            pairs+=1
    rows=[{"n":n,"principal":principal[n],"nonprincipal":other[n],
           "trivial_character":even[n],
           "nontrivial_character":odd[n]} for n in range(1,min(limit,30)+1)]
    return {
      "status":"exact finite integer arithmetic; known ideal-class identification not re-proved",
      "field":"Q(sqrt(-5))", "discriminant":-20,
      "forms":["x^2+5y^2","2x^2+2xy+3y^2"],
      "class_group":"C2", "up_to_n":limit,
      "coprime_multiplicativity_pairs_tested":pairs,
      "nonmultiplicative_principal_witness":{"f(2)":0,"f(3)":0,"f(6)":2},
      "character_projection":"a_+ = principal+nonprincipal, a_- = principal-nonprincipal",
      "first_coefficients":rows,
      "scope":"finite character algebra and counting; Hecke Euler products are prior art; no RH"
    }

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--limit",type=int,default=128)
    a=p.parse_args()
    if a.limit>10000:
        p.error("limit <=10000")
    print(json.dumps(report(a.limit),indent=2))
