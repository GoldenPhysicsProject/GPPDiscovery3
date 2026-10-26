"""A positive prime Hodge complex with EXACT signed Weil prime supertrace.

For each p set q=p**(-1/2), d_p=S-qI between two copies of l2(N0),
with S the unilateral isometry. Its range is closed, ker d_p=0, and
H^1=coker d_p is one-dimensional, spanned by the normalized thermal
vector v_q=(sqrt(1-q²)q^n)_{n>=0}. The shift S commutes with d_p and
acts on H^1 by the positive number q; in odd cohomological degree its
SUPERTRACE is -q^k for S^k. Thus the full signed von Mangoldt prime
distribution is a supertrace of positive-metric, genuine cohomology.

No RH proof: the positive Archimedean completion and an arithmetic
Hodge/OS even-cohomology theorem are still required.
"""

import math

import numpy as np

from prime_os_weil_boundary import autocorrelation, primes_through


def finite_cokernel_test(p: int, n: int = 60, max_k: int = 6):
    q=p**(-0.5)
    d=np.zeros((n+1,n))
    for j in range(n):
        d[j,j]=-q
        d[j+1,j]=1.0
    v=np.array([q**j for j in range(n+1)])
    v/=np.linalg.norm(v)
    # v spans coker(d), which is exactly 1D for this rectangular section.
    assert np.linalg.norm(d.T@v)<1e-12
    s=np.diag(np.ones(n),k=-1)
    result=[]
    sk=np.eye(n+1)
    for k in range(1,max_k+1):
        sk=s@sk
        measured=float(v@sk@v)
        expected=q**k
        exact_finite=q**k*(1-q**(2*(n+1-k)))/(1-q**(2*(n+1)))
        assert abs(measured-exact_finite)<1e-12
        assert abs(measured-expected)<1e-8
        result.append((k,measured,-math.log(p)*measured))
    return result


def signed_arithmetic_supertrace(h: float):
    """Exact all-prime *cohomological* supertrace on compact F=g*gtilde."""
    total=0.0
    terms=0
    for p in primes_through(math.ceil(math.exp(2*h))):
        q=p**(-0.5)
        n=p
        k=1
        while math.log(n)<2*h:
            f=autocorrelation(math.log(n),h)
            if f:
                total-=math.log(p)*f*q**k
                terms+=1
            n*=p
            k+=1
    return {"half_support":h,"signed_cohomology_supertrace":total,
            "weil_prime_quadratic_contribution":2*total,
            "prime_power_terms":terms}


def archimedean_graded_density(x: float):
    """v34 real-place density, as an exact vacuum-minus-odd-tower trace.

    For x>0, w_inf(x)=exp(x/2)-sum_{j>=0}exp(-(5/2+2j)x).
    After half-density conjugation the tower has POSITIVE frequencies
    3,5,7,... and an even vacuum at frequency zero.
    """
    assert x>0
    original=math.exp(-x/2)+math.exp(x/2)-math.exp(-x/2)/(-math.expm1(-2*x))
    graded=math.exp(x/2)-math.exp(-2.5*x)/(-math.expm1(-2*x))
    conjugated=1-math.exp(-3*x)/(-math.expm1(-2*x))
    assert math.isclose(original,graded,rel_tol=1e-11,abs_tol=1e-11)
    assert math.isclose(math.exp(-x/2)*original,conjugated,
                        rel_tol=1e-11,abs_tol=1e-11)
    return {"x":x,"w_infty":original,
            "half_density_graded_vacuum_minus_gamma":conjugated}


if __name__=="__main__":
    for p in (2,3,5,11):
        values=finite_cokernel_test(p,n=60)
        print(f"p={p} k=1 measured={values[0][1]:.14g} "
              f"supertrace_logp={values[0][2]:.14g} "
              f"k=6 measured={values[-1][1]:.14g}")
    for h in (1.0,1.5,2.0,2.5,3.0):
        d=signed_arithmetic_supertrace(h)
        print(" ".join(f"{key}={value:.14g}" if isinstance(value,float)
                       else f"{key}={value}" for key,value in d.items()))
    for x in (0.1,0.3,math.log(2),1.0,2.0):
        d=archimedean_graded_density(x)
        print("arch "+" ".join(f"{key}={value:.14g}" for key,value in d.items()))
