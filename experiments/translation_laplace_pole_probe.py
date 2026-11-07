#!/usr/bin/env python3
"""Floating control for the fixed-profile *continuous-translation* RH gate.

Exact underlying theorem (not established by this script): choose a real
compactly supported bump train V with entire Laplace W(z)=integral Ve^{zu}
and no zeros in Re(z)>0. Translate V to V_a(u)=V(u-a). The Weil sewing
scalar is C(a)=sum_{z=rho-1/2}m_rho W(z)^2 exp(2*a*z).

RH iff C(a) is bounded on a>=a0: under RH every exponential has modulus
one; conversely, if bounded its Laplace transform in a is holomorphic on
Re(p)>0. For Re(p)>1 it equals sum m_rho W(z)^2 exp(-(p-2z)a0)/(p-2z),
which extends meromorphically and has a genuine pole at p=2z for EVERY
off-line zero with Re(z)>0 (W(z)!=0). Hence boundedness excludes them all,
without any rightmost-zero or phase-dominance assumption.

The numerical simulation uses only the first N known on-line zeta zeros
and, separately, the deliberately shifted off-critical F_1/4 control.
It does not establish boundedness of the actual full zeta arithmetic form.
"""

import argparse
import cmath
import json
import math


def first_zero_modes(k: int, theta: float, n_zeros: int):
    import mpmath as mp
    mp.mp.dps = 42
    delta=mp.mpf(7)/10
    roots=[mp.exp(-delta*(2*m+mp.mpf('.5'))) for m in range(k-1)]
    def profile(z):
        s=2*mp.pi**2*mp.expm1(z)/(z*(z*z+4*mp.pi**2))
        return s*mp.fprod(mp.exp(delta*z)-r for r in roots)
    modes=[]
    for n in range(1,n_zeros+1):
        gamma=mp.im(mp.zetazero(n))
        for shift in ([mp.mpf(0)] if theta==0 else
                      [-mp.mpf(1)/4,mp.mpf(1)/4]):
            z=shift+1j*gamma
            modes.append((complex(z),complex(profile(z)**2)))
    return modes


def signal(a: float,modes) -> float:
    return 2*math.fsum((c*cmath.exp(2*a*z)).real for z,c in modes)


def rms_window(center: float,modes) -> float:
    points=(signal(center-.5+j/300,modes) for j in range(301))
    return math.sqrt(math.fsum(x*x for x in points)/301)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--k',type=int,default=6)
    ap.add_argument('--zeros',type=int,default=60)
    args=ap.parse_args()
    if args.k<2 or args.k>10 or args.zeros<10 or args.zeros>150:
        ap.error('K=2..10 and zero pairs=10..150')
    centers=[.35,5.,10.,15.,20.,25.]
    results=[]
    for theta in (0.,.25):
        modes=first_zero_modes(args.k,theta,args.zeros)
        rms=[rms_window(c,modes) for c in centers]
        if theta==0:
            assert max(rms)/min(rms)<2
        else:
            # This is a finite-zero planted-control stress test only.
            assert rms[-1]>100*rms[0]
        results.append({'theta':theta,'K':args.k,'zero_pairs':args.zeros,
                        'centers':centers,'sewing_rms_numeric':rms})
    print(json.dumps({'status':'Numerical finite-zero translation control; '
                   'analytic boundedness equivalence is RH-strength and OPEN',
                   'results':results},indent=2))


if __name__=='__main__':
    main()
