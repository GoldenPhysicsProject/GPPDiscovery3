"""Exact linear-field OS/Gram completion of the signed Weil prime source.

The Gram field is linear in g, so its norm is quadratic in g.  Only the
prime powers whose shifts overlap a compactly supported test enter; the
contact term is the exact cost of writing the negative prime off-diagonal
as a sum of squares.  This is a candidate bulk-to-boundary accounting
identity, not the whole Weil Archimedean completion.
"""

import math

from prime_os_weil_boundary import autocorrelation, primes_through


def linear_prime_bonds(h: float) -> dict[str, float | int]:
    """For support g subset [-h,h], calculate sum w ||g-U_l g||^2.

    w=Lambda(p**k)/sqrt(p**k), l=log(p**k).  Translating g by l defines
    a linear map to a direct sum of L2(R) copies, so its norm square is
    OS-positive independently of any zero information.
    """
    norm2=autocorrelation(0.0,h)
    b=0.0
    negative_prime=0.0
    terms=0
    for p in primes_through(math.ceil(math.exp(2*h))):
        n=p
        while math.log(n)<2*h:
            logn=math.log(n)
            w=math.log(p)/math.sqrt(n)
            b+=w
            negative_prime-=2*w*autocorrelation(logn,h)
            terms+=1
            n*=p
    contact=2*b*norm2
    norm_square=contact+negative_prime
    return {"test_half_support":h,"prime_power_terms":terms,
            "g_norm_square":norm2,"contact_required":contact,
            "signed_weil_prime":negative_prime,
            "prime_bond_norm_square":norm_square,
            "contact_to_prime_abs":contact/abs(negative_prime)}


if __name__=="__main__":
    for h in (1,1.5,2,2.5,3):
        result=linear_prime_bonds(h)
        assert result["prime_bond_norm_square"] >= -1e-10
        print(" ".join(f"{k}={v:.12g}" if isinstance(v,float) else f"{k}={v}"
                       for k,v in result.items()))
