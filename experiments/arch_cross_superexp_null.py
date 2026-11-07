#!/usr/bin/env python3
"""Constructive all-K near-null witness for the sine² cross-Arch Hankel Gram.

THEOREM (analytic, not merely the printed finite checks). Let K>=2,
a_j=7/20+(7/10)j, lambda_n=2n+1/2, and
  M_K[j,k]=sum_(n>=0) L_n(j)L_n(k),
  L_n(j)=exp(-lambda_n*a_j)*S(-lambda_n),
  S(-lambda)=integral_0^1 sin(pi*u)^2 exp(-lambda*u) du.
Then M_K is strictly positive definite, but

  0 < lambda_min(M_K) <= exp(-X_K)/(Lambda_K^2*(1-exp(-7/5))),
  Lambda_K=2K-3/2,
  X_K=(7/5)*(K-1)*(K-3/2)+(7/10)*Lambda_K.

PROOF. Set r_n=exp(-(7/10)*lambda_n) and choose the nonzero vector
of coefficients of P_K(z)=prod_(m=0..K-2)(z-r_m). Its top coefficient
is 1, so ||v||^2>=1. Since P_K(r_n)=0 for n<K-1 and
|P_K(r_n)|<=prod_(m<K-1) r_m for n>=K-1, the Rayleigh quotient of
M_K is <= (prod r_m)^2 * sum_(n>=K-1) exp(-2a_0*lambda_n)/lambda_n^2.
Use 0<S(-lambda)<=1/lambda and sum the geometric tail. This gives
the stated bound. Strict positive definiteness follows because a
polynomial of degree <=K-1 cannot vanish at K distinct r_n.

Numerical displays are non-rigorous; the inequalities tagged EXACT use
only rational Taylor lower bounds exp(X)>sum_(j<=D) X^j/j! and
exp(7/5)>1+7/5+(7/5)^2/2>2.

No claim about full prime/pole/diagonal Arch completion or RH.
"""

import json
from fractions import Fraction as Q
from math import factorial, log10


def exact_upper(k: int, degree: int = 140) -> Q:
    if not 2 <= k <= 50:
        raise ValueError("check K>=2; executable bound limited to K<=50")
    a0, delta = Q(7, 20), Q(7, 10)
    lam = Q(2*k, 1) - Q(3, 2)
    x = 2*delta*(k-1)*(Q(k)-Q(3, 2)) + 2*a0*lam
    exp_lower = sum((x**j / factorial(j) for j in range(degree+1)), Q(0))
    assert 1 + Q(7, 5) + Q(7, 5)**2/2 > 2
    return Q(2) / (lam**2 * exp_lower)


def main() -> None:
    k6 = exact_upper(6)
    k10 = exact_upper(10)
    assert k6 < Q(1, 10**18)
    assert k10 < Q(1, 10**52)
    # The decimal approximations are for reading only; comparisons above
    # are done with unrounded rational integers.
    print(json.dumps({
        "theorem": "For every K>=2: 0<lambda_min(-C_cross,K)<=exp(-X_K)/(Lambda_K^2*(1-exp(-7/5)))",
        "source": "Wave234 folded sine2 CROSS-arch only",
        "test_K6": {"certified_less_than": "1e-18", "upper_approx": f"{float(k6):.7e}"},
        "test_K10": {"certified_less_than": "1e-52", "upper_approx": f"{float(k10):.7e}"},
        "consequence": "Strict finite-K definiteness coexists with a super-exponentially disappearing minimal spectral gap.",
        "rh": "not established"
    }, indent=2))


if __name__ == "__main__":
    main()
