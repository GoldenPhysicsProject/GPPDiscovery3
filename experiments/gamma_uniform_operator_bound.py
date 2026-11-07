#!/usr/bin/env python3
"""All-K operator-norm bound for the PURE Arch gamma Toeplitz + Hankel blocks.

For unit-length sine2 bumps at centres a_j=a0+delta*j, delta=7/10,
the gamma Toeplitz entries T_ij=A_|i-j| obey, for d>=2,

  A_d=-sum_lambda exp[-lambda*(delta*d-1)] S(-lambda)^2,
  |A_d|<=|A_2| exp[-delta*lambda_min*(d-2)].

Consequently uniformly over ALL sizes K,

  ||T_K||_op <= |A0|+2|A1|+2|A2|/(1-exp(-delta*lambda_min)).

The folded cross-Arch Hankel block is a negative Gram matrix

  H_ij=-sum_lambda h_lambda(i)h_lambda(j),
  h_lambda(j)=exp[-lambda*(a0+delta*j)] S(-lambda),

so ||H_K||_op <= Tr(-H_K) <=
  sum_lambda exp(-2*lambda*a0) S(-lambda)^2 /
             (1-exp(-2*lambda*delta)) < infinity.

Thus ||T_K+H_K|| is bounded independent of K.  The atomic prime sewing
C-channel in Turbo Wave306 has linearly growing sections; PURE gamma
Arch alone cannot absorb those.  This does not exclude compensation
from prime/pole/SQ or non-additive completed sewing, and does not prove
RH.  The operator inequality is exact; printed decimals use floating
compact quadrature for A0,A1 and are NOT certified enclosures.

The cross-trace infinite remainder after N modes has an optional exact
rational cap via S(-lambda)<20/lambda^3, elementary exponential Taylor
lower bounds, and geometric tails.
"""

import argparse
from fractions import Fraction
import json
import math
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent))
from arch_annihilator_full_s1_probe import A0, DELTA, sine2_laplace
from gamma_offsupport_laplace import (
    compact_gamma_toeplitz, laplace_offsupport, rational_tail_cap,
    rational_taylor_exp,
)


def cross_trace_tail_cap(theta: Fraction, n_modes: int = 80,
                         degree: int = 36) -> Fraction:
    """Outward RATIONAL upper bound for infinite Gram trace after N."""
    if theta not in (Fraction(0),Fraction(1,4)) or n_modes<1:
        raise ValueError("theta 0 or 1/4; n_modes positive")
    a0,delta=Fraction(7,20),Fraction(7,10)
    offsets=([Fraction(1,2)] if theta==0 else
             [Fraction(1,2)-theta,Fraction(1,2)+theta])
    qtail_upper=(1+4*a0)/(4*a0)
    result=Fraction(0)
    for off in offsets:
        lam=2*n_modes+off
        geometric_upper=(1+2*lam*delta)/(2*lam*delta)
        exp_upper=1/rational_taylor_exp(2*lam*a0,degree)
        result+=Fraction(400)*exp_upper*qtail_upper*geometric_upper/lam**6
    return result


def gamma_constants(theta: float, n_modes: int = 100) -> dict:
    """Only finite-gamma constants are numerical; K-independence is exact."""
    if theta not in (0., .25):
        raise ValueError("theta=0 or 1/4")
    shifts=(0.,) if theta==0 else (-theta,theta)
    lam_min=.5-theta
    q=math.exp(-DELTA*lam_min)
    a0=compact_gamma_toeplitz(0,theta)
    a1=compact_gamma_toeplitz(1,theta)
    a2=laplace_offsupport(2,theta,n_modes)
    exact_theta=Fraction(0) if theta==0 else Fraction(1,4)
    a2_tail=rational_tail_cap(2,exact_theta,n_modes)
    # These are approximate decimals because a0/a1 require certified
    # finite quadrature.  The omitted d=2 infinite tail IS enclosed.
    toeplitz_bound=abs(a0)+2*abs(a1)+2*(abs(a2)+float(a2_tail))/(1-q)
    trace=math.fsum(
        math.exp(-2*lam*A0)*sine2_laplace(-lam)**2/
        (1-math.exp(-2*DELTA*lam))
        for sh in shifts for n in range(n_modes)
        for lam in (2*n+.5+sh,))
    cap=cross_trace_tail_cap(exact_theta,n_modes)
    return {"theta":theta,"gamma_lags":(a0,a1,a2),"lambda_min":lam_min,
            "decay_ratio_q":q,"toeplitz_norm_upper_numeric":toeplitz_bound,
            "a2_mode_tail_exact_upper_fraction":str(a2_tail),
            "a2_mode_tail_upper_float":float(a2_tail),
            "cross_trace_head_numeric":trace,
            "cross_trace_tail_exact_upper_fraction":str(cap),
            "cross_trace_tail_upper_float":float(cap),
            "combined_uniform_norm_bound_numeric":toeplitz_bound+trace+float(cap)}


def matrix_check(k: int, theta: float, n_modes: int = 80) -> dict:
    """Optional numpy check at finite K, never an all-K certification."""
    import numpy as np
    if k<2 or k>128:
        raise ValueError("2<=K<=128")
    vals=[compact_gamma_toeplitz(d,theta) if d<=1 else
          laplace_offsupport(d,theta,n_modes) for d in range(k)]
    t=np.array([[vals[abs(i-j)] for j in range(k)] for i in range(k)])
    h=np.zeros((k,k))
    js=np.arange(k)
    for sh in ((0.,) if theta==0 else (-theta,theta)):
        for n in range(n_modes):
            lam=2*n+.5+sh
            v=np.exp(-lam*(A0+DELTA*js))*sine2_laplace(-lam)
            h-=np.outer(v,v)
    total=t+h
    return {"K":k,"theta":theta,
            "norm_toeplitz":float(np.linalg.norm(t,ord=2)),
            "norm_hankel":float(np.linalg.norm(h,ord=2)),
            "norm_combined":float(np.linalg.norm(total,ord=2)),
            "hankel_largest_eigenvalue":float(np.linalg.eigvalsh(h)[-1])}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n-modes',type=int,default=80)
    p.add_argument('--matrix-sizes',type=int,nargs='*',default=[6,16,40])
    args=p.parse_args()
    results=[]
    for theta in (0.,.25):
        row=gamma_constants(theta,args.n_modes)
        row['finite_matrix_checks']=[matrix_check(k,theta,args.n_modes)
                                     for k in args.matrix_sizes]
        for x in row['finite_matrix_checks']:
            assert x['norm_combined'] < row['combined_uniform_norm_bound_numeric']+1e-9
        results.append(row)
    print(json.dumps({
        "status":"PAPER all-K pure Arch operator bound; numerical constants and finite matrix checks, not RH",
        "results":results},indent=2))


if __name__=='__main__':
    main()
