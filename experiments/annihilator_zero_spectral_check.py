#!/usr/bin/env python3
"""Independent numerical prime/zero explicit-formula bridge for Wave234 S1.

For compactly supported V_K, put F_K(z)=integral V_K(u) exp(z*u) du.
The exact Weil explicit formula gives on zeta's two-sector operator
Q(V_K)=sum_{nontrivial zeta zeros rho} F_K(rho-1/2)F_K(1/2-rho),
C(V_K)=sum_{nontrivial zeta zeros rho} F_K(rho-1/2)^2.
These identities are unconditional; the numeric experiment samples only the
first N *known on-critical-line* conjugate zero pairs and is NOT a proof that
all zeros are on the line, nor a rigorous infinite-tail enclosure.

For F_theta= zeta(s-theta)zeta(s+theta), replace rho by rho +/- theta;
these known shifted zeros form a controlled off-line-zero adversary.

Requires mpmath>=1.3 and scipy/numpy (optional lab check only).
"""

import argparse
import json


def zero_amplitudes(k: int, theta: float, ordinates: list) -> tuple[float,float]:
    import mpmath as mp
    a0, delta = mp.mpf(7)/20, mp.mpf(7)/10
    roots = [mp.exp(-delta*(2*m+mp.mpf('.5'))) for m in range(k-1)]

    def bump_mellin(z):
        S = 2*mp.pi**2*mp.expm1(z)/(z*(z*z+4*mp.pi**2))
        return (mp.exp(a0*z)*S
                *mp.fprod(mp.exp(delta*z)-root for root in roots))

    shifts = [mp.mpf(0)] if theta==0 else [-mp.mpf(1)/4,mp.mpf(1)/4]
    q,c = mp.mpf(0),mp.mpf(0)
    for y in ordinates:
        for shift in shifts:
            z=shift+1j*y
            a=bump_mellin(z)
            b=bump_mellin(-z)
            q+=2*mp.re(a*b)
            c+=2*mp.re(a*a)
    return float(q),float(c)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--zeros',type=int,default=120)
    ap.add_argument('--cutoff',type=float,default=500.)
    args=ap.parse_args()
    if not 10<=args.zeros<=400:
        ap.error('restrict to 10..400 for bounded CI runtime')

    import mpmath as mp
    from arch_annihilator_full_s1_probe import evaluate
    mp.mp.dps=45
    ordinates=[mp.im(mp.zetazero(n)) for n in range(1,args.zeros+1)]
    cases=[]
    for k in (6,10):
        for theta in (0.,.25):
            zero_q,zero_c=zero_amplitudes(k,theta,ordinates)
            source=evaluate(k,theta,args.cutoff)
            n2=source['norm2']
            pos=source['parity']['plus']['full_over_norm']
            neg=source['parity']['minus']['full_over_norm']
            prime_q=.5*n2*(pos+neg)
            prime_c=.5*n2*(pos-neg)
            delta_q=abs(prime_q-zero_q)
            delta_c=abs(prime_c-zero_c)
            if args.zeros>=120:
                # Empirical consistency test, NOT a certified zero-tail bound.
                assert delta_q<1e-8, (k,theta,'diagonal',delta_q)
                assert delta_c<1e-8, (k,theta,'sewing',delta_c)
            cases.append({
                'K':k,'theta':theta,'zeros_pairs':args.zeros,
                'gamma_t_cutoff':args.cutoff,
                'prime_pole_arch_diag_numeric':prime_q,
                'first_N_zeros_diag_numeric':zero_q,
                'diagonal_residual':delta_q,
                'prime_pole_arch_sewing_numeric':prime_c,
                'first_N_zeros_sewing_numeric':zero_c,
                'sewing_residual':delta_c})
    print(json.dumps({
        'status':'Numerical cross-check of unconditional explicit-formula identities;'
                 ' finite known-line zeros, not RH verification',
        'cases':cases},indent=2))


if __name__=='__main__':
    main()
