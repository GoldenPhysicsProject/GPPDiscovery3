"""Character-flavor controls on the thermal prime-field OS construction.

For chi_4 the arithmetic coefficients vary in sign, but the connected
prime-field Gram covariance remains a positive innovation square. This
demonstrates that negative arithmetic contributions can be physical
interference and that bulk OS positivity alone is insufficient for GRH.
"""

from __future__ import annotations

import math

from prime_os_weil_boundary import autocorrelation, primes_through


def chi4(p: int) -> int:
    return 0 if p==2 else (1 if p%4==1 else -1)


def arithmetic_character_field(h: float) -> dict[str,float|int]:
    means={"trivial":0.0,"chi4":0.0}
    variances={"trivial":0.0,"chi4":0.0}
    count=0
    for p in primes_through(math.ceil(math.exp(2*h))):
        data=[]
        n=p
        k=1
        while math.log(n)<2*h:
            a=math.log(p)*autocorrelation(math.log(n),h)
            if a:
                data.append((k,a))
                count+=1
            n*=p
            k+=1
        q=p**(-0.5)
        for name,character in (("trivial",1),("chi4",chi4(p))):
            coeffs={k:a*character**k for k,a in data}
            means[name]+=sum(a*q**k for k,a in coeffs.items())
            variances[name]+=(1-q*q)*sum(
                abs(sum(a*q**(k-j) for k,a in coeffs.items() if k>=j))**2
                for j in range(1,max(coeffs,default=0)+1))
    return {"h":h,"prime_power_terms":count,
            "zeta_source":means["trivial"],
            "chi4_source":means["chi4"],
            "zeta_connected_variance":variances["trivial"],
            "chi4_connected_variance":variances["chi4"]}


if __name__=="__main__":
    for h in (1,1.5,2,2.5,3):
        d=arithmetic_character_field(h)
        assert d["zeta_connected_variance"]>=0
        assert d["chi4_connected_variance"]>=0
        print(" ".join(f"{key}={value:.11g}" if isinstance(value,float)
                       else f"{key}={value}" for key,value in d.items()))
