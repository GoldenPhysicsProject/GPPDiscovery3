#!/usr/bin/env python3
"""Exact point-count-integrality falsifier for shifted zeta F_theta.

Only compares LOCAL Euler denominators with what ordinary finite-field
point counts permit. Does not prove RH or geometric realization of zeta.
"""
from fractions import Fraction as Q
import json

def exact_gate():
    # At p=2, theta=1/4:
    # Z_fake(T) = 1/[(1-2^(1/4)T)(1-2^(-1/4)T)].
    # Its coefficient of T must be #X(F_2) if from a variety over F_2.
    # Write x=2^(1/4), so 1<x<5/4 since (5/4)^4>2.
    # AM-GM, strictly, gives x+1/x>2; while x+1/x<x+1<9/4<3.
    upper=Q(5,4)+1
    assert Q(5,4)**4>2
    assert Q(2)<upper<Q(3)
    # The coefficient is therefore strictly between consecutive integers.
    # In contrast, P1(F2) has 2+1=3 points; for every p,
    # Z(P1/Fp,T)=1/[(1-T)(1-pT)] has integral first coefficient p+1.
    return {
        "prime":2,"theta":"1/4",
        "fake_first_point_count":"2^(1/4)+2^(-1/4)",
        "proved_strict_lower":"2", "proved_strict_upper":"9/4",
        "point_count_integral":False,
        "sphere_first_point_count":3,
        "scope":"Ftheta is not itself an ordinary finite-field zeta with these Euler factors; neither RH nor BSD follows"
    }
if __name__=="__main__":
    print(json.dumps(exact_gate(),sort_keys=True,indent=2))
