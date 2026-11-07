#!/usr/bin/env python3
"""Directed Arb enclosures of the actual off-support gamma Toeplitz entries.

For x=7d/10 >= 1 (integer d>=2), the digamma/Fourier identity proved
in gamma_offsupport_laplace.py gives the unconditional exact real series

  A_gamma(x) = -sum_(lambda) exp[-lambda*(x-1)] S(-lambda)^2,
  S(-lambda)=2*pi^2*(1-exp(-lambda))/(lambda*(lambda^2+4*pi^2)),
  lambda=2n+1/2 (zeta), or lambda=2n+1/2 +/- 1/4 (Ftheta).

Every finite term is enclosed by rigorous 384-bit Arb operations. The
entire omitted series is bounded by the EXACT Fraction inequality in
gamma_offsupport_laplace.rational_tail_cap. Adding symmetric uncertainty
of the rational upper cap gives a certified interval, serialized as
EXACT dyadic endpoints (integer mantissa, binary exponent), not decimal
strings whose displayed digits might conceal rounding.

This closes gamma off-support entry enclosure under the analytic source
identity. The compact d=0,1 quadratures and the rest of the arithmetic
matrix remain separate. No RH claim.

Requires python-flint==0.9.0; default catalog run is disabled and the
optional Arb CI lane installs this dependency.
"""

import argparse
from fractions import Fraction
import json
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent))
from gamma_offsupport_laplace import rational_tail_cap


def certificate(d: int, theta: Fraction, n_modes: int = 60,
                taylor_degree: int = 36, precision_bits: int = 384) -> dict:
    from flint import arb, ctx, fmpq

    if d<2 or theta not in (Fraction(0),Fraction(1,4)):
        raise ValueError("d>=2 and theta 0 or 1/4 required")
    if n_modes<1 or precision_bits<128:
        raise ValueError("n_modes>=1 and precision>=128 required")
    ctx.prec=precision_bits
    ar=lambda q:arb(fmpq(q.numerator,q.denominator))

    L=ar(Fraction(7*d,10)-1)
    shifts=([Fraction(0)] if theta==0 else [-theta,theta])
    pi2=arb.pi()**2
    partial=arb(0)
    for shift in shifts:
        for n in range(n_modes):
            lam=ar(Fraction(4*n+1,2)+shift)
            s=2*pi2*(1-(-lam).exp())/(lam*(lam*lam+4*pi2))
            partial-=(-lam*L).exp()*s*s

    cap=rational_tail_cap(d,theta,n_modes,taylor_degree)
    # Arb tracks roundoff in 'partial'. arb(0,error.upper()) encloses
    # *every* possible tail with |tail|<=cap without discarding that
    # tracked roundoff. Its radius is rounded safely upwards.
    error_upper=ar(cap).upper()
    enclosed=partial+arb(0,error_upper)
    assert enclosed < 0, (d,theta,enclosed)

    lower_m,lower_e=enclosed.lower().man_exp()
    upper_m,upper_e=enclosed.upper().man_exp()
    return {
        "d":d,"theta":str(theta),"x":str(Fraction(7*d,10)),
        "modes_per_gamma_shift":n_modes,"precision_bits":precision_bits,
        "tail_fraction":str(cap),"tail_upper_float_for_display":float(cap),
        "certified_negative":True,
        "dyadic_lower":{"mantissa":str(lower_m),"exponent":int(lower_e)},
        "dyadic_upper":{"mantissa":str(upper_m),"exponent":int(upper_e)},
        "ball_display_only":enclosed.str(35),
    }


def exact_dyadic(point: dict) -> Fraction:
    """Recover an endpoint with exact integer arithmetic, no float."""
    return Fraction(int(point['mantissa']))*Fraction(2)**int(point['exponent'])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--d-max',type=int,default=5)
    p.add_argument('--n-modes',type=int,default=60)
    p.add_argument('--precision-bits',type=int,default=384)
    p.add_argument('--cross-precision-bits',type=int,default=512,
                   help="independent precision rerun; 0 disables overlap audit")
    p.add_argument('--taylor-degree',type=int,default=36)
    args=p.parse_args()
    if not 2<=args.d_max<=100:
        p.error("2<=d-max<=100")
    rows=[]
    for theta in (Fraction(0),Fraction(1,4)):
        for d in range(2,args.d_max+1):
            row=certificate(d,theta,args.n_modes,args.taylor_degree,
                            args.precision_bits)
            if args.cross_precision_bits:
                again=certificate(d,theta,args.n_modes,args.taylor_degree,
                                  args.cross_precision_bits)
                lo1=exact_dyadic(row['dyadic_lower'])
                hi1=exact_dyadic(row['dyadic_upper'])
                lo2=exact_dyadic(again['dyadic_lower'])
                hi2=exact_dyadic(again['dyadic_upper'])
                assert lo1<=hi2 and lo2<=hi1, ('disjoint enclosures',d,theta)
                row['verified_cross_precision_bits']=args.cross_precision_bits
                row['certified_overlap']=True
            rows.append(row)
    print(json.dumps({
        "status":"CERTIFIED gamma off-support Laplace series under paper-derived source identity; no RH claim",
        "interval_representation":"EXACT signed dyadic [mantissa*2^exponent] endpoints",
        "rows":rows},indent=2))


if __name__=='__main__':
    main()
