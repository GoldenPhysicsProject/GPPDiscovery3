#!/usr/bin/env python3
"""Celestial CP1 nodal intersections vs false principal-circle intersection.

stdlib geometry and arithmetic tests; --numeric requires optional mpmath.
No proof RH or that celestial CFT Hilbert space equals arithmetic spectrum.
"""
import argparse
import cmath
import json
import math


def cayley(s):
    """s -> w=(s-1)/s; s=0 is the pole."""
    return (s-1)/s


def latitude(w):
    """Inverse stereographic celestial n_z coordinate."""
    return (1-abs(w)**2)/(1+abs(w)**2)


def dimension_circle(d,s):
    """Delta_d=d*s: w_d=(Delta_d-d)/Delta_d equals w_1(s)."""
    if d<=0:raise ValueError("d>0")
    return (d*s-d)/(d*s)


def prime_phase(p,t):
    return cmath.exp(-1j*t*math.log(p))


def eta_partial(s,n):
    return sum(((-1)**(k-1))*cmath.exp(-s*math.log(k))
               for k in range(1,n+1))


def numerical_probe(count=4,dps=45):
    """Independent theta completion integral and first critical-zero nodes."""
    import mpmath as mp
    mp.mp.dps=dps
    def H(u):
        x=mp.exp(u)
        return mp.fsum(mp.exp(-mp.pi*n*n*x) for n in range(1,7))
    def theta_integral(t):
        return mp.quad(lambda u:H(u)*mp.exp(u/4)*mp.cos(t*u/2),
                       [0,mp.mpf('.5'),1,2,4,10])
    def xi(s):
        return mp.mpf('.5')*s*(s-1)*mp.power(mp.pi,-s/2)*mp.gamma(s/2)*mp.zeta(s)
    zeros=[]
    theta=mp.mpf('0.25')
    for k in range(1,count+1):
        rho=mp.zetazero(k);t=mp.im(rho);w=(rho-1)/rho
        I=theta_integral(t);wave=mp.mpf('.5')-(t*t+mp.mpf('.25'))*I
        shifted=[]
        for sign in (-1,1):
            s=rho+sign*theta
            ww=(s-1)/s
            lat=(1-abs(ww)**2)/(1+abs(ww)**2)
            shifted.append({'Re_s':mp.nstr(mp.re(s),8),
                'n_z':mp.nstr(lat,23),
                'abs_Ftheta':mp.nstr(abs(xi(s-theta)*xi(s+theta)),9)})
        zeros.append({'index':k,'t':mp.nstr(t,27),'w':mp.nstr(w,24),
            'abs_xi':mp.nstr(abs(xi(rho)),8),
            'theta_integral':mp.nstr(I,27),
            'target_integral':mp.nstr(1/(2*(t*t+mp.mpf('.25'))),27),
            'theta_wave_residual':mp.nstr(abs(wave),8),
            'fake_shifted_off_equator':shifted})
    meet=2*mp.pi/mp.log(mp.mpf(3)/2)
    return {'eta_half_at_t0':mp.nstr(mp.altzeta(mp.mpf('.5')),25),
        'prime_2_3_phase_meet_t':mp.nstr(meet,25),
        'zeta_at_meet':mp.nstr(mp.zeta(mp.mpf('.5')+1j*meet),20),
        'interpretation':'numerical replay only; genuine RH and dynamical intertwiner open',
        'zeros':zeros}


def report():
    s=.5+14.13472514173469j
    w=cayley(s)
    return {'w':str(w),'equator_latitude':latitude(w),
        'same_circle_d_1_to_6':
         all(abs(dimension_circle(d,s)-w)<1e-12 for d in range(1,7)),
        'interpretation':'all dimension-normalized circles coincide at ALL parameters'}


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--numeric',action='store_true')
    ap.add_argument('--count',type=int,default=4)
    ap.add_argument('--dps',type=int,default=45)
    a=ap.parse_args()
    print(json.dumps(numerical_probe(a.count,a.dps) if a.numeric else report(),indent=2))
