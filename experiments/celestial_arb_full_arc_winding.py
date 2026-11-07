#!/usr/bin/env python3
"""Potential certified argument-principle count on a height-centered celestial circle.

Uses python-flint==0.9.0 Arb complex balls, not floating samples.
ALL contour ARCS must be enclosed and EXCLUDE ZERO before winding is trusted.
If any segment enclosure fails, this script exits NONZERO; no certification.
This proves only a finite local count, NEVER RH or global spectral positivity.
"""
import argparse
import json


def verify(segments=2048,bits=256):
    from flint import acb, arb, ctx
    if segments<128 or bits<128:
        raise ValueError("segments>=128, bits>=128")
    ctx.prec=bits
    a=arb(1)/2
    r=arb(9)/20
    shift=arb(1)/4
    h=arb("14.1347251417346937904572519835624702707842571156992431756855")
    pi=arb.pi()
    I=acb(0,1)
    # fixed rational center h is distinct from the separately Arb-certified
    # height gamma. The width of the first zero interval must be controlled.
    rho=acb.zeta_zeros(1,1)[0]
    gamma=rho.imag
    if not abs(gamma-h)<arb("1e-30"):
        raise ArithmeticError("first certified Arb critical zero not close enough to fixed rational h", gamma,h)
    if not rho.real==arb(1)/2:
        raise ArithmeticError("first zero enclosure not on the classical certified line")

    def s_of_theta(theta):
        w=r*(I*theta).exp()
        return acb(arb(1)/2)+I*h+a*(1+w)/(1-w)

    def xi(s):
        s=acb(s)
        return (s*(s-1)/2) * acb.pi()**(-s/2) * (s/2).gamma() * s.zeta()

    sampled_true=[]
    sampled_fake=[]
    min_arc_true=None
    min_arc_fake=None
    for k in range(segments):
        center=(arb(2*k+1)*pi)/segments
        halfwidth=pi/segments
        t=arb(center,halfwidth) # rigorously covers ENTIRE angular subarc
        s=s_of_theta(t)
        xi_arc=xi(s)
        fake_arc=xi(s-shift)*xi(s+shift)
        if xi_arc.contains(0):
            raise ArithmeticError("cannot certify xi zero-free ARC",k,segments,str(xi_arc))
        if fake_arc.contains(0):
            raise ArithmeticError("cannot certify fake zero-free ARC",k,segments,str(fake_arc))
        lo1=xi_arc.abs_lower()
        lo2=fake_arc.abs_lower()
        if min_arc_true is None or lo1<min_arc_true:min_arc_true=lo1
        if min_arc_fake is None or lo2<min_arc_fake:min_arc_fake=lo2

        tpoint=(arb(2*k)*pi)/segments
        sp=s_of_theta(tpoint)
        sampled_true.append(xi(sp))
        sampled_fake.append(xi(sp-shift)*xi(sp+shift))

    def winding(vals,expected,label):
        total=arb(0)
        maxangle=arb(0)
        for k in range(segments):
            ratio=vals[(k+1)%segments]/vals[k]
            angle=ratio.arg()
            if not (angle<pi/2 and angle>-pi/2):
                raise ArithmeticError("phase increment not certified within (-pi/2,pi/2)",label,k,str(angle))
            total+=angle
            if abs(angle)>maxangle:maxangle=abs(angle)
        turns=total/(2*pi)
        if not (turns.contains(expected) and
                turns>arb(expected)-arb(1)/4 and
                turns<arb(expected)+arb(1)/4):
            raise ArithmeticError("winding integer not uniquely certified",label,str(turns),expected)
        return {"certified_winding":expected,
                "turns_interval":turns.str(30),
                "max_endpoint_phase_increment":maxangle.str(15)}

    z=winding(sampled_true,0,"xi")
    f=winding(sampled_fake,1,"F1/4")
    return {
        "status":"RIGOROUS_ARB_COMPLETE_ARC_WINDING_CERTIFICATE",
        "segments":segments,"precision_bits":bits,
        "chart_height_fixed_rational":h.str(42),"chart_scale":"1/2",
        "contour_radius":"9/20",
        "first_known_critical_zero_height_enclosure":gamma.str(42),
        "xi":z,"Ftheta":f,
        "lower_abs_arc_xi":min_arc_true.str(20),
        "lower_abs_arc_Ftheta":min_arc_fake.str(20),
        "scope":"Only one specified celestial disk, not RH or Weil positivity",
        "method":"Each full analytic image arc enclosed by a convex acb rectangle excluding zero; winding of endpoint polygon is homotopic to actual curve; Arb argument increments bound a unique integer."
    }


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--segments",type=int,default=2048)
    ap.add_argument("--bits",type=int,default=256)
    x=ap.parse_args()
    print(json.dumps(verify(x.segments,x.bits),indent=2))
