#!/usr/bin/env python3
"""Exact rational uniform tail estimate for the folded sine-squared arch sewing.

Setup in the GPP-bridge Wave234 S1 sector: bumps sin(pi*(u-a_j))**2
on [a_j,a_j+1], a_j=7/20+7*j/10. Its arch cross-kernel is
  -exp(-v/2)/(1-exp(-2*v)) = -sum_{n>=0}exp(-(2*n+1/2)*v).
Consequently C_arch[j,k]=-sum_{n>=0}L_n(j)*L_n(k), where
  L_n(j)=exp(-lam_n*a_j)*S(-lam_n), lam_n=2*n+1/2,
  S(-lam)=2*pi**2*(1-exp(-lam))/(lam*(lam**2+4*pi**2)).

For n>=8, |S(-lam_n)|<=20/lam_n**3, using pi**2<10.
The common tail is bounded by
  400/(16**6) * exp(-231/20)/(1-exp(-7/5))
  < 800/(16**6*50000) < 1e-9.

Only rational finite Taylor inequalities are checked here. The input
identity for the arch kernel and the elementary inequality pi**2<10 are
part of the paper proof, not computed from floating approximations.
This applies ONLY to the stated cross-arch channel, not to the full
arithmetic 6x6 matrix or the all-K RH margin.
"""

from fractions import Fraction
from math import factorial


def exp_lower(x: Fraction, degree: int) -> Fraction:
    assert x >= 0
    return sum((x**k / factorial(k) for k in range(degree + 1)), Fraction(0))


def main() -> None:
    # Both statements have proofs using ONLY exact rational arithmetic;
    # exp(x) exceeds every positive Taylor partial sum for x>0.
    t1 = exp_lower(Fraction(231, 20), 12)
    t2 = exp_lower(Fraction(7, 5), 2)
    assert t1 > 50_000, "need exp(231/20)>50000"
    assert t2 > 2, "need exp(7/5)>2"

    tail_bound = Fraction(800, 16**6 * 50_000)
    assert tail_bound < Fraction(1, 10**9)
    print("exp(231/20) > Taylor_12 > 50000: EXACT PASS")
    print("exp(7/5) > Taylor_2 > 2: EXACT PASS")
    print(f"uniform cross-arch tail for every j,k, after n=0..7: "
          f"< {tail_bound} < 1/10^9")
    print("No statement about diagonal-arch, prime, pole, or full matrix error")


if __name__ == "__main__":
    main()
