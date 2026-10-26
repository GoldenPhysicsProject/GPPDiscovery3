"""Exact thermal-circle OS construction for two opposite prime-field orientations.

Prime p is an independent bosonic oscillator species, E_p=log(p). For beta=1
the local KMS oscillator has thermal Schwinger 2-point function
  S_p(tau) = [p^(-tau)+p^(-(1-tau))]/(1-p^(-1)), 0 <= tau <= 1.
Reflection about tau=0 yields a positive two-feature Gram on 0<t<1/2.
Higher free Schwinger functions follow from Wick pairing. The infinite product
state exists on the quasi-local algebra via consistent finite-prime states.

No claim that this Gaussian OS system equals Weil's signed arithmetic form.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np


PRIMES=(2,3,5,7,11,13,17,19,23,29)


def thermal_covariance(p, tau):
    if not 0 <= tau <= 1:
        raise ValueError("thermal separation must be in [0,1]")
    return (p**(-tau)+p**(tau-1))/(1-p**(-1))


def feature(p, t):
    if not 0 <= t <= 0.5:
        raise ValueError("reflection half-circle t in [0,1/2]")
    return np.array([p**(-t),p**(t-0.5)])/math.sqrt(1-1/p)


def thermal_os_matrix(p, times):
    a=np.array([feature(p,t) for t in times])
    return a@a.T


def wick_four(S2):
    """Centered free prime-field fourth moment from pairings."""
    return S2(0,1)*S2(2,3)+S2(0,2)*S2(1,3)+S2(0,3)*S2(1,2)


def tfd_shift_expectation_exact(p, k):
    """Prove coefficient algebraically via geometric sum:

    <Omega| (S+^dagger S-^dagger)^k |Omega>
      =(1-q^2) sum_{n>=0} q^(2n+k)=q^k, q=1/sqrt(p).
    """
    return p ** (-k/2)


def finite_tfd_shift_cutoff(p,k,N):
    """Unnormalized partial geometric expectation; converges to q^k."""
    q=p**(-.5)
    return (1-q*q)*sum(q**(2*n+k) for n in range(N))


def main():
    times=np.array([.035,.12,.21,.32,.48])
    print("Prime thermal OS reflection: kernel rank <=2 and PSD")
    gram=np.zeros((len(times),len(times)))
    for p in PRIMES:
        gp=thermal_os_matrix(p,times)
        for i,t in enumerate(times):
            for j,u in enumerate(times):
                assert abs(gp[i,j]-thermal_covariance(p,t+u))<1e-12
        assert min(np.linalg.eigvalsh(gp))>-1e-10
        gram+=gp
        print(f"p={p:2d} S_p(.3)={thermal_covariance(p,.3):.8f} "
              f"reflection-eig-max={max(np.linalg.eigvalsh(gp)):.8f}")
    print("all 10 prime species combined min eigenvalue",min(np.linalg.eigvalsh(gram)))

    # Wick-normal-ordered n-particle OS Gram matrices are n! times the
    # entrywise nth power of the one-field reflected covariance. The Schur
    # product theorem makes them PSD for all n; test degrees 1..5 here.
    print("Wick/Schwinger hierarchy: normal-ordered n-particle OS positivity")
    for degree in range(1,6):
        gram_n=math.factorial(degree)*np.power(gram,degree)
        smallest=min(np.linalg.eigvalsh(gram_n))
        assert smallest > -1e-8
        print(f"Wick degree={degree}, minimum OS eigenvalue={smallest:.9g}")

    print("TFD prime-power oriented-shift vacuum amplitudes")
    for p in (2,3,5):
        for k in (1,2,3):
            got=finite_tfd_shift_cutoff(p,k,300)
            want=tfd_shift_expectation_exact(p,k)
            assert abs(got-want)<1e-13
            print(f"p={p} k={k} amplitude={want:.12f}; desired Lambda(p^k)/sqrt(p^k) "
                  f"=log(p)*amplitude={math.log(p)*want:.12f}")

    print("Exact local thermal occupation, from normalized geometric distribution")
    for p in PRIMES[:5]:
        mean=1/(p-1)
        fluct=(p)/(p-1)**2
        assert mean>0 and fluct>0
        print(f"p={p} mean={mean:.9f} variance={fluct:.9f}")

    # Exact prime-power von Mangoldt midpoint identity. It directly relates
    # the arithmetic log-derivative source to the two-orientation KMS
    # Schwinger function on the critical thermal circle.
    print("Thermal midpoint reproduces every local von Mangoldt prime-power sum")
    for p in PRIMES:
        q=p**(-.5)
        lhs=(math.log(p)/2)*(1+q)*thermal_covariance(p,.5)
        rhs=math.log(p)*q/(1-q)
        assert abs(lhs-rhs)<1e-12
        print(f"p={p:2d}: midpoint source={lhs:.12f}, Euler power sum={rhs:.12f}")


if __name__ == "__main__":
    main()
