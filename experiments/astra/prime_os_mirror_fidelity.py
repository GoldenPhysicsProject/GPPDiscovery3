"""All-prime mirror-KMS vacuum disjointness: a sharp critical-line selection.

For 0<sigma<1, the local normalized doubled oscillator vector is
Omega_{sigma,t,p}=sqrt(1-p^{-2sigma}) sum_{k>=0}
   p^{-sigma*k} exp(-it*k*log p)|k,k>.
Reflect by the *anti-Hermitian* critical real structure
rho=sigma+it -> 1-conj(rho)=1-sigma+it, keeping spectral phase t.

The local Hellinger amplitude is exactly
R_p=sqrt((1-p^{-2sigma})(1-p^{-2(1-sigma)}))/(1-p^{-1}).
It is 1 at sigma=1/2; for sigma != 1/2 its all-prime product is ZERO.
The result uses only divergence of sum_p 1/p, so no zeta-zero assumption.

Consequence: any bounded quasi-local operator has vanishing infinite-
volume matrix elements between distinct mirror KMS sectors, so the
two sectors cannot be coupled by local prime-field interference.
This DOES NOT prove RH: it is still necessary to prove that each
zeta zero must arise from a nonzero physical matrix element between
the mirror sectors, with the full Archimedean scattering/completion.
"""

from __future__ import annotations

import math

from prime_os_weil_boundary import primes_through


def local_overlap(p: int, sigma: float) -> float:
    assert 0<sigma<1
    a=p**(-sigma)
    b=p**(-(1-sigma))
    r=math.sqrt((1-a*a)*(1-b*b))/(1-a*b)
    assert math.isclose(r*r,1-((a-b)/(1-a*b))**2,
                        rel_tol=1e-12,abs_tol=1e-12)
    return r


def fidelity_to_prime(limit:int,sigma:float):
    ps=primes_through(limit)
    log_product=sum(math.log(local_overlap(p,sigma)) for p in ps)
    # Quantified universal proof bound (valid all p):
    # log R_p<=-0.5*((a-b)/(1-1/p))**2.
    proof_bound=-0.5*sum(
        ((p**(-sigma)-p**(-(1-sigma)))/(1-1/p))**2 for p in ps)
    assert log_product<=proof_bound+1e-9
    return {"sigma":sigma,"prime_limit":limit,"prime_species":len(ps),
            "vacuum_fidelity":math.exp(log_product),
            "log_fidelity":log_product,"upper_bound_log":proof_bound}


def mirror_fidelity_susceptibility(limit:int):
    """Exact -d²/dtheta² log fidelity at theta=sigma-1/2=0.

    chi_P=4 sum_{p<=P} (log p)**2/[p(1-1/p)**2]. By PNT this
    grows ~2(log P)**2. The coherent width is of order 1/log P.
    """
    return 4*sum(math.log(p)**2/(p*(1-1/p)**2)
                 for p in primes_through(limit))


if __name__=="__main__":
    for sigma in (0.3,0.4,0.45,0.5,0.55,0.6,0.7):
        for cutoff in (100,1000,10000,100000):
            r=fidelity_to_prime(cutoff,sigma)
            print(f"sigma={sigma:.2f} primes<={cutoff}: "
                  f"fidelity={r['vacuum_fidelity']:.8g}, "
                  f"log={r['log_fidelity']:.7g}, bound={r['upper_bound_log']:.7g}")
    for cutoff in (100,1000,10000,100000):
        chi=mirror_fidelity_susceptibility(cutoff)
        r=fidelity_to_prime(cutoff,0.49)
        gaussian=-0.5*chi*(0.01**2)
        print(f"P={cutoff} susceptibility={chi:.9g} "
              f"logF(theta=.01)={r['log_fidelity']:.9g} "
              f"quadratic_prediction={gaussian:.9g}")
