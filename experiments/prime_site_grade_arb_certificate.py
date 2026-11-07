#!/usr/bin/env python3
"""Directed Arb certificate for a p=3, K=3 prime-site signed parity form.

Inherits ONLY the explicit Wave234 S1 sine² G0 and local prime-power source
identity. All transcendental operations (log, exp, cos, sin, pi) are enclosed
at 384 and independently at 512 bits. All prime-power support is bounded
EXACTLY by G0's |argument|<1 condition. Certificates are exact dyadic
intervals. Proves determinant/minor SIGNS of this one finite source matrix,
not positivity of the global prime-Arch form, not any RH theorem.

Requires python-flint==0.9.0 (optional Arb GitHub CI lane).
"""
from fractions import Fraction as Q
import json


def _dyad_interval(a):
    lm, le = a.lower().man_exp()
    um, ue = a.upper().man_exp()
    return (Q(int(lm)) * Q(2)**int(le),
            Q(int(um)) * Q(2)**int(ue))


def _dyad_str(a):
    lm, le = a.lower().man_exp()
    um, ue = a.upper().man_exp()
    return {"lower": {"mantissa": str(lm), "exponent": int(le)},
            "upper": {"mantissa": str(um), "exponent": int(ue)}}


def certificate(bits=384):
    from flint import arb, ctx, fmpq
    ctx.prec=bits
    def aq(a):
        a=Q(a)
        return arb(fmpq(a.numerator,a.denominator))

    zero=arb(0)
    one=arb(1)
    pi=arb.pi()
    logp=arb(3).log()
    delta=aq(Q(7,10)); a0=aq(Q(7,20))
    maxlog=aq(Q(11,2))   # 2(a0+delta*(K-1))+2 with K=3

    # The source's G0 is even, supported in [-1,1].
    def G0(x):
        u=abs(x)
        if u>=one:
            return zero
        if not (u<one):
            raise ArithmeticError("support boundary unresolved by Arb",u)
        phase=2*pi*u
        return ((one-u)*(2+phase.cos())/8
                +3*phase.sin()/(16*pi))

    terms=[]
    for a in range(1,7):
        t=a*logp
        if t<maxlog:
            terms.append((t, logp*(-t/2).exp()))
        elif not t>maxlog:
            raise ArithmeticError("source support cutoff unresolved",a,t)
    if len(terms)!=5:
        raise ArithmeticError("unexpected number of p=3 powers",len(terms))

    def build(ftheta):
        def kernel(h):
            total=zero
            for t,w in terms:
                if ftheta:
                    factor=(t/4).exp()+(-t/4).exp()
                else:
                    factor=one
                total+=w*factor*G0(t-h)
            return total
        d=[kernel(delta*j)+kernel(-delta*j) for j in range(3)]
        c=[kernel(2*a0+delta*l+1) for l in range(5)]
        # sigma=-1: -(D-C).
        M=[[-(d[abs(i-j)]-c[i+j]) for j in range(3)] for i in range(3)]
        a,b,c0=M[0]
        d0,e,f=M[1]
        g,h,i0=M[2]
        minors=[a, a*e-b*d0,
                a*(e*i0-f*h)-b*(d0*i0-f*g)+c0*(d0*h-e*g)]
        if not (minors[0]>zero and minors[1]>zero):
            raise ArithmeticError("first two minors not certified >0",minors)
        if ftheta:
            if not minors[2]>zero:
                raise ArithmeticError("Ftheta determinant not certified >0",minors[2])
        else:
            if not minors[2]<zero:
                raise ArithmeticError("zeta determinant not certified <0",minors[2])
        return M,minors

    zM,zD=build(False)
    fM,fD=build(True)
    return {
        "bits":bits,
        "source":"K3 p3 sigma-minus Wave234 S1 prime-only, full p site",
        "zeta_inertia":"(2 positive, 1 negative)",
        "ftheta_inertia":"(3 positive, 0 negative)",
        "zeta_minor_intervals":[_dyad_str(x) for x in zD],
        "ftheta_minor_intervals":[_dyad_str(x) for x in fD],
        "internal_matrix_balls":[zM,fM],
        "internal_minor_balls":[zD,fD],
    }


def run():
    lo=certificate(384)
    hi=certificate(512)
    # Assert independently computed 384/512-bit source *matrix entries*
    # and leading minors overlap, with exact rational dyadic comparisons.
    def overlaps(a,b):
        a0,a1=_dyad_interval(a)
        b0,b1=_dyad_interval(b)
        return a0<=b1 and b0<=a1
    for key in ("internal_matrix_balls", "internal_minor_balls"):
        xs=lo[key]; ys=hi[key]
        for A,B in zip(xs,ys):
            if key=="internal_matrix_balls":
                pairs=[(A[i][j],B[i][j]) for i in range(3) for j in range(3)]
            else:
                pairs=list(zip(A,B))
            if not all(overlaps(a,b) for a,b in pairs):
                raise AssertionError("independent Arb interval overlap failed",key)
    return {
        "status":"DIRECTED_ARBITRARY_PRECISION_INTERVAL_CERTIFIED",
        "rho_zeros_assumed":False,"RH_claimed":False,
        "support_powers":"p^1..p^5, source G0 compact support",
        "inertia_zeta":lo["zeta_inertia"],"inertia_ftheta":lo["ftheta_inertia"],
        "precision_bits":[384,512],
        "zeta_minors_384":lo["zeta_minor_intervals"],
        "ftheta_minors_384":lo["ftheta_minor_intervals"],
        "zeta_minors_512":hi["zeta_minor_intervals"],
        "ftheta_minors_512":hi["ftheta_minor_intervals"],
        "warning":"finite site matrix only: no poles/gamma, no global positivity"
    }


if __name__=="__main__":
    print(json.dumps(run(),indent=2))
