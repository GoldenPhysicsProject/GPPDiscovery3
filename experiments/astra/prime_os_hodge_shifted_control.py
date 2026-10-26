"""Hostile shifted-zeta control for the positive prime Hodge source.

For |theta|<1/2, a pair of positive local Hodge cohomologies with
q_{p,+}=p^{-1/2+theta} and q_{p,-}=p^{-1/2-theta} has odd supertrace
-2 cosh(theta*k*log p) Lambda(p**k)/sqrt(p**k).

This is the signed prime term for zeta(s-theta)zeta(s+theta), which has
zeros off the usual critical line whenever theta != 0 (since zeta has
critical zeros). Thus prime-local positivity, even with orientation
pairing, CANNOT be the final RH selection rule; the completed
Archimedean/boundary dynamics must break this countermodel.
"""

from __future__ import annotations

import math

from prime_os_weil_boundary import autocorrelation, primes_through


def shifted_pair(h: float, theta: float):
    assert abs(theta)<0.5
    total_hodge=0.0
    direct=0.0
    smallest_gap=1.0
    count=0
    for p in primes_through(math.ceil(math.exp(2*h))):
        qp=p**(-0.5+theta)
        qm=p**(-0.5-theta)
        assert 0<qp<1 and 0<qm<1
        smallest_gap=min(smallest_gap,1-max(qp,qm))
        n=p
        k=1
        while math.log(n)<2*h:
            f=autocorrelation(math.log(n),h)
            if f:
                # The Hodge negative sign appears only in the supertrace.
                total_hodge-=math.log(p)*(qp**k+qm**k)*f
                direct-=2*math.log(p)/math.sqrt(n)*math.cosh(
                    theta*math.log(n))*f
                count+=1
            n*=p
            k+=1
    assert math.isclose(total_hodge,direct,rel_tol=1e-13,abs_tol=1e-13)
    return {"h":h,"theta":theta,"terms":count,
            "odd_positive_cohomology_supertrace":total_hodge,
            "independently_evaluated_shifted_prime_source":direct,
            "uniform_local_hodge_gap_lower_bound":smallest_gap}


if __name__=="__main__":
    for theta in (0.0,0.1,0.2,0.4):
        d=shifted_pair(h=3.0,theta=theta)
        print(" ".join(f"{key}={value:.12g}" if isinstance(value,float)
                       else f"{key}={value}" for key,value in d.items()))
