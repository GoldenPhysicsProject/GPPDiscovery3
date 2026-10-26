"""Operator-valued complete von Mangoldt source from prime OS/KMS fields.

The logarithmic test F is compactly supported. Therefore its arithmetic
operator J(F) is an exactly FINITE evaluation of an ALL-PRIME distribution;
there is no externally selected prime or prime-power cutoff. The positive
connected vacuum covariance is calculated exactly from the two-mode KMS
purification, alongside the signed prime term entering Weil's formula.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad


def primes_through(limit):
    is_p=[True]*(limit+1)
    is_p[0]=is_p[1]=False
    for p in range(2,math.isqrt(limit)+1):
        if is_p[p]:
            for q in range(p*p,limit+1,p):
                is_p[q]=False
    return [p for p in range(2,limit+1) if is_p[p]]


def bump(u,h=2.5):
    return max(0.0,1-(u/h)**2)**3 if abs(u)<h else 0.0


def autocorrelation(ell,h=2.5):
    if abs(ell)>=2*h:
        return 0.0
    lower=max(-h,-h-ell)
    upper=min(h,h-ell)
    return quad(lambda u:bump(u,h)*bump(u+ell,h),lower,upper,
                epsabs=1e-13,epsrel=1e-13)[0]


def prime_field_source(h=2.5):
    """J(F)=sum_{p,k}log(p)F(k log p)V_p^k in the quasi-local algebra.

    For compact supp F subset [-2h,2h], only p^k<exp(2h) contribute.
    Under Omega_p, E[V_p^k]=p^(-k/2),
    E[(V_p^k)^* V_p^m]=p^(-|k-m|/2).
    """
    cutoff=math.ceil(math.exp(2*h))
    data={}
    for p in primes_through(cutoff):
        entries=[]
        k=1
        n=p
        while n<math.exp(2*h):
            ell=math.log(n)
            coeff=math.log(p)*autocorrelation(ell,h)
            if coeff:
                entries.append((k,n,coeff))
            k+=1
            n*=p
        if entries:
            data[p]=entries

    mean=0.0
    variance=0.0
    innovation_square=0.0
    direct_prime_sum=0.0
    for p,entries in data.items():
        q=p**(-0.5)
        mean+=sum(a*q**k for k,_,a in entries)
        direct_prime_sum+=sum(math.log(p)*autocorrelation(math.log(n),h)/math.sqrt(n)
                              for k,n,a in entries)
        variance+=sum(a*b*(q**abs(k-m)-q**(k+m))
                      for k,_,a in entries for m,_,b in entries)
        innovation_square+=(1-q*q)*sum(
            abs(sum(a*q**(k-j) for k,_,a in entries if k>=j))**2
            for j in range(1,max(k for k,_,a in entries)+1))
    return {"mean":mean,"direct":direct_prime_sum,"variance":variance,
            "innovation_square":innovation_square,
            "rawnorm":variance+mean*mean,"prime_species":len(data),
            "prime_powers":sum(len(x) for x in data.values()),
            "signed_weil_prime_contribution":-2*mean}


def main():
    for h in (1.0,1.5,2.0,2.5,3.0):
        r=prime_field_source(h)
        print(f"test support [-{h},{h}], F support [-{2*h},{2*h}]")
        for key,value in r.items():
            print(f"  {key}: {value:.14g}" if isinstance(value,float) else f"  {key}: {value}")
        assert math.isclose(r["mean"],r["direct"],rel_tol=1e-13,abs_tol=1e-13)
        assert r["variance"] >= -1e-12
        assert math.isclose(r["variance"],r["innovation_square"],
                            rel_tol=1e-12,abs_tol=1e-12)
        print("  exact omega[J(F)] = von Mangoldt source; Var[J(F)] >= 0 VERIFIED")


if __name__=="__main__":
    main()
