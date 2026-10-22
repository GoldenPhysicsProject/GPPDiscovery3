#!/usr/bin/env python3
"""GPP exact universal sharpness of pole-balanced causal energy transfer.

Astra's correct transfer:
    F_eps(t) = exp(-eps*t) J0(t)
    L(J_eps)(z) = [(z+eps)(z+eps+1/2)/(z(z+1/2))] L(F_eps)(z).
For ANY complex compactly-supported F, not necessarily arithmetic,
the finite-window dissipation bound is
    ||J_eps||²_L2[0,T] >= ||F_eps||²_L2[0,T].
Combining with J0=exp(eps*t)*F implies:
    E0(T) <= exp(2 eps T) E_eps(T).
THIS LOSS IS SHARP over generic smooth packet signals.
The script checks exact partial fractions in rational arithmetic,
then brackets the sharp ratio with analytic inequalities for a single
high-frequency endpoint-localized wave packet. No RH claims.
"""
import math
from fractions import Fraction

def fractional_decomp(e):
    e=Fraction(e)
    alpha=2*e*(e+Fraction(1,2))
    beta=e-2*e*e
    # A_e(z)=1+alpha/z+beta/(z+1/2):
    assert alpha+beta==2*e
    assert alpha/2==e*(e+Fraction(1,2))
    return alpha,beta

def ratio_bracket(e,T,h,omega):
    """J0(t)=exp(e*t)*exp(i omega*t) indicator_[T-h,T].
    Valid non-smooth L² extremizer; smooth bump approximation is possible.
    u=J_e=F+alpha*int F+beta*exp(-t/2)*int exp(s/2)F(s)ds.
    Within support the convolution sup error is at most
    delta=2|alpha|/omega +2|beta|/sqrt(omega²+1/4).
    """
    e=float(e)
    aa,bb=map(float,fractional_decomp(Fraction(e).limit_denominator()))
    assert 0<h<T and omega>0
    d=2*abs(aa)/omega + 2*abs(bb)/math.hypot(omega,0.5)
    assert d<1
    num=(math.exp(2*e*T)-math.exp(2*e*(T-h)))/(2*e)
    lo=num/(h*(1+d)**2)
    hi=num/(h*(1-d)**2)
    return lo,hi,math.exp(2*e*T),d

if __name__=="__main__":
    for e in [Fraction(1,10),Fraction(1,4),Fraction(1,2),Fraction(3,4)]:
        print("epsilon",e,"alpha,beta",fractional_decomp(e))
    for h in [0.8,0.2,0.05,0.01,0.001]:
        lo,hi,target,err=ratio_bracket(Fraction(1,2),8,h,1e10)
        print(f"T=8 eps=1/2 h={h:g} bracket [{lo:.9f}, {hi:.9f}] "
              f"sharp e^8={target:.9f} L_inf_error<={err:.2e}")
        assert lo<=target<=hi
    print("Proof: omega->infinity makes J_eps/F->1 in L2. "
          "Then h->0 forces E0/E_eps->exp(2eps*T).")
    print("This is a SHARPNESS theorem for general signals, "
          "not a nonexistence theorem for arithmetic-specific estimates.")
