#!/usr/bin/env python3
"""GPP 2026-10-23: critical prime TFD two-time clock scaling (NO RH).

In the finite-prime critical paired TFD state
 Omega_{P,1}=tensor_{p<=P} sqrt(1-p^-1) sum_{k>=0} p^-k/2 |k,k>,
the sum boost H_+ + H_- has energy 2 sum_{p<=P} k_p log p.

At fixed physical t != 0, the overlap with its sum-time evolution
tends to zero as P tends to infinity (classical zero-free line,
no RH assumption); the fixed-t state is in a disjoint ITP sector.

But X_P = (H_+ + H_-)/(2 log P) converges in law to the
Dickman distribution with characteristic exponent
 int_0^1 (e^(it u)-1)/u du, and cumulants kappa_m=1/m.
The alternate Gibbs beta=1+delta global normalized state has
delta log n converging to Exp(1), with cumulants (m-1)!.

THIS IS NOT a zeta-zero detector. Both limits are prime-pole/PNT
phenomena and do not assert anything about nontrivial zeros.
"""
import math


def primes_upto(limit):
    composite=bytearray(limit+1)
    for n in range(2,int(math.isqrt(limit))+1):
        if not composite[n]:
            for j in range(n*n,limit+1,n):
                composite[j]=1
    return (p for p in range(2,limit+1) if not composite[p])


def cumulants_prime_cutoff(P):
    L=math.log(P)
    c1=c2=c3=c4=0.0
    for p in primes_upto(P):
        w=math.log(p)/L
        q=1.0/p
        v=1-q
        c1+=w*q/v
        c2+=w*w*q/v**2
        c3+=w**3*q*(1+q)/v**3
        c4+=w**4*q*(1+4*q+q*q)/v**4
    return (c1,c2,c3,c4)


if __name__=="__main__":
    print("P       k1        k2        k3        k4")
    for P in (1000,10000,100000,1000000):
        vals=cumulants_prime_cutoff(P)
        print(f"{P:<7}",*(f"{v:10.7f}" for v in vals))
    print("Dickman limit:",*[1/m for m in range(1,5)])
    print("Gibbs pole/Exp(1) limit:",*[math.factorial(m-1) for m in range(1,5)])
    assert abs(cumulants_prime_cutoff(1000000)[2]-1/3)<1e-3
    assert abs(cumulants_prime_cutoff(1000000)[1]-1/2)<2e-3
    print("PASS: finite prime-cutoff third and second cumulants approach Dickman limits")
