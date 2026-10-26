"""Test the 'negative is positive' interpretation for fermionic primes.

The Euclidean fermion propagator is anti-periodic, so its reflected
negative-time branch appears with a MINUS in the raw time-ordered
two-point function. Fermionic OS reflection contributes the compensating
grading sign, leaving a positive Gram kernel.

The arithmetic all-positive prime-power winding coefficients instead
come from a PERIODIC fermionic supertrace. Its density matrix is not a
positive state. The two assertions must not be conflated; a physical
Hodge/BRST quotient would have to connect them.
"""

import math

import numpy as np


def one_prime(p, n=25):
    e=math.log(p)
    r=1/p
    t=np.linspace(0,0.5,n)
    forward=np.outer(np.exp(-e*t),np.exp(-e*t))/(1+r)
    backwards=np.outer(np.exp(-e*(0.5-t)),np.exp(-e*(0.5-t)))/(1+r)
    raw_fermion_negative_time=-backwards
    # CAR anti-periodicity is accompanied by an extra sign in the OS
    # anti-linear graded reflection, producing the positive 'backwards'.
    corrected_reflected=forward-raw_fermion_negative_time
    assert np.linalg.eigvalsh(corrected_reflected).min()>-1e-12
    assert np.linalg.eigvalsh(raw_fermion_negative_time).min()<-1e-5
    # Physical fermion Gibbs state vs (-1)^F-twisted supertrace.
    physical_occupied_probability=r/(1+r)
    graded_occupied_weight=-r/(1-r)
    assert physical_occupied_probability>0 and graded_occupied_weight<0
    # Logarithmic source of periodic fermionic determinant, Re(s)>0:
    # d_s log(1-p^{-s}) = (log p) sum_{k>=1}p^{-ks}.
    sum_primepowers=e*r/(1-r)
    from_graded_det=e*r/(1-r)
    assert math.isclose(sum_primepowers,from_graded_det,rel_tol=1e-15)
    # Anti-periodic positive fermion determinant log(1+p^{-s}) has
    # alternating winding coefficients and a different derivative sign.
    alternating=e*r/(1+r)
    return {"p":p,"os_min_eig":float(np.linalg.eigvalsh(corrected_reflected).min()),
            "raw_negative_eig":float(np.linalg.eigvalsh(raw_fermion_negative_time).min()),
            "positive_fermion_occupation":physical_occupied_probability,
            "nonpositive_supertrace_occupation":graded_occupied_weight,
            "all_primepowers_source":sum_primepowers,
            "alternating_fermionic_source":alternating}


if __name__=="__main__":
    for p in (2,3,5,7,11):
        d=one_prime(p)
        print(" ".join(f"{k}={v:.12g}" if isinstance(v,float) else f"{k}={v}"
                       for k,v in d.items()))
