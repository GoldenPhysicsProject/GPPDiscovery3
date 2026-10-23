#!/usr/bin/env python3
"""2026-10-21 GPP: critical GCD/Haar Möbius Euler packet cannot control prefix energy.

Completely exact Fraction tests (NO zeta zeros, NO RH), prime sets with
first 1..7 primes.  Infinite proof:
 F_P=prod_{p in P}(1-chi_p/sqrt(p)) in H_GCD (critical Poisson).
 ||F_P||^2=prod_{p in P}(1-1/p) -->0.
 On [0,log 2), prefix cumulative amplitude is 1, so the target
 half-line energy ∫ e^{-t}|sum_{d|primorial,d<=e^t} mu(d)/sqrt(d)|²dt >=1/2.
 Hence the map taking coefficient polynomials to this arithmetic
 prefix-energy space has no bounded extension from H_GCD.
 Unbounded does NOT imply unclosable; RH not proved.

 At weights d^{-(1/2+theta)}, the critical-Haar local norm is
 1-2/p^{1+theta}+1/p^{1+2theta}; product tends
 0 at theta=0, positive finite at theta>0, infinity at -1/2<theta<0.
"""
from fractions import Fraction
from itertools import combinations
from math import gcd
from functools import reduce
from operator import mul

PRIMES=(2,3,5,7,11,13,17)

def divisors_with_mu(primes):
    out=[(1,1)]
    for p in primes:
        out=out+[(d*p,-sgn) for d,sgn in out]
    return out

def lcm(m,n):
    return m//gcd(m,n)*n

def gcd_norm_exact(primes):
    terms=divisors_with_mu(primes)
    return sum((Fraction(sm*sn,lcm(m,n)) for m,sm in terms for n,sn in terms),Fraction(0))

def euler_product_exact(primes):
    return reduce(mul,(Fraction(p-1,p) for p in primes),Fraction(1))

def check():
    for k in range(1,len(PRIMES)+1):
        ps=PRIMES[:k]
        lhs=gcd_norm_exact(ps)
        rhs=euler_product_exact(ps)
        assert lhs==rhs,(ps,lhs,rhs)
        target_vacuum=Fraction(1,2)
        ratio_lower_bound=target_vacuum/lhs
        print(f"{k} primes, P={reduce(mul,ps)}: GCD norm²={lhs} "
              f"≈{float(lhs):.10g}; target energy≥1/2; ratio≥{float(ratio_lower_bound):.9g}")
    print("PASS: exact rational identity for first 1..7 primes.")
    print("Infinite theorem uses Euler divergence sum_p 1/p and the fixed first cell [1,2).")
    print("Unboundedness of natural transfer is proved; closability is NOT addressed.")
if __name__=="__main__":
    check()
