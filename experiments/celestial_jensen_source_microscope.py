#!/usr/bin/env python3
"""Local celestial Cayley microscope: Jensen defect and argument-principle winding.

EXACT: half-plane/disk map; shifted-control interior-root coordinate;
Jensen lower bound conditional on known zeta zero; positivity of theta
Fourier kernel and its shifted convolution (documented separately).
NUMERICAL: finite contour winding / quadrature, NOT zero certification.
Optional mpmath required only for --numeric. Does not prove RH.
"""
import argparse
import cmath
import json
import math


def local_cayley(s, t0, a):
    """Map Re(s)>1/2 onto unit disk, recentering at height t0."""
    if a <= 0:
        raise ValueError("scale must be >0")
    z = s - (0.5 + 1j*t0)
    return (z-a)/(z+a)


def local_inverse(w, t0, a):
    if a <= 0 or abs(w) >= 1:
        raise ValueError("a>0, |w|<1 needed for local inverse")
    return .5+1j*t0+a*(1+w)/(1-w)


def local_latitude(s,t0,a):
    w=local_cayley(s,t0,a)
    return (1-abs(w)**2)/(1+abs(w)**2)


def fake_root_coordinate(theta,a):
    """For shifted zero s=1/2+theta+i*t0, the local disk coordinate."""
    if not 0<theta<a:
        raise ValueError("0<theta<a needed")
    return (theta-a)/(theta+a)


def jensen_lower_bound(r,root_radius):
    if not 0<root_radius<r<1:
        raise ValueError("0<|root|<r<1")
    return math.log(r/root_radius)


def winding_samples(values):
    """Phase-lift along a sampled *closed* contour; no certification."""
    if len(values)<4 or abs(values[0]-values[-1])>1e-7*max(1,abs(values[0])):
        raise ValueError("closed contour samples required")
    if any(z==0 for z in values):
        raise ValueError("zero exactly on sampled contour")
    increments=[cmath.phase(values[k+1]/values[k])
                for k in range(len(values)-1)]
    if max(abs(v) for v in increments)>=math.pi*.75:
        raise ValueError("phase undersampling; refine contour")
    return sum(increments)/(2*math.pi),max(abs(v) for v in increments)


def jensen_samples(values,center):
    """Mean log modulus; equals exact Jensen defect only in continuum."""
    if center==0 or any(z==0 for z in values):
        raise ValueError("nonzero center/contour required")
    return sum(math.log(abs(v/center)) for v in values)/len(values)


def numerical_probe(dps=38,samples=160):
    import mpmath as mp
    if dps<25 or samples<48:
        raise ValueError("dps>=25, samples>=48")
    mp.mp.dps=dps
    gamma=mp.im(mp.zetazero(1))
    a=mp.mpf('.5'); theta=mp.mpf('.25'); r=mp.mpf('.45')
    def s_of_w(w):
        return mp.mpf('.5')+1j*gamma+a*(1+w)/(1-w)
    def xi(s):
        return s*(s-1)/2*mp.power(mp.pi,-s/2)*mp.gamma(s/2)*mp.zeta(s)
    def fake(s):
        return xi(s-theta)*xi(s+theta)
    def sample(fun):
        v0=fun(s_of_w(0))
        vv=[fun(s_of_w(r*mp.exp(2*mp.pi*1j*j/samples)))
            for j in range(samples)]
        # mp.arg handles exact phase to high working precision, then 53-bit floats
        wraps=[float(mp.arg(vv[(j+1)%samples]/vv[j]))
               for j in range(samples)]
        if max(abs(v) for v in wraps)>=2.35:
            raise AssertionError("contour undersampled")
        return {"approx_winding":sum(wraps)/(2*math.pi),
                "approx_jensen":float(mp.fsum(mp.log(abs(v/v0)) for v in vv)/samples),
                "min_contour_modulus":float(min(abs(v) for v in vv)),
                "max_phase_increment":max(abs(v) for v in wraps)}
    def H(u):
        return mp.fsum(mp.exp(-mp.pi*n*n*mp.exp(u)) for n in range(1,8))
    def theta_source(s):
        I=mp.quad(lambda u: H(u)*mp.exp(u/4)*
                  mp.cosh((s-mp.mpf('.5'))*u/2),
                  [0,mp.mpf('.5'),1,2,4,10])
        return mp.mpf('.5')+s*(s-1)*I
    def phi_pos(x):
        x=abs(x)
        return mp.fsum((4*mp.pi**2*n**4*mp.exp(mp.mpf('4.5')*x)-
                        6*mp.pi*n*n*mp.exp(mp.mpf('2.5')*x))*
                       mp.exp(-mp.pi*n*n*mp.exp(2*x)) for n in range(1,7))
    # Previous project theta-tilt convolution; positive for every real u.
    def conv_positive(u):
        return mp.quad(lambda x:phi_pos(x)*phi_pos(u-x)*
                       mp.cosh(theta*(u-2*x)),[-4,-2,-1,0,1,2,4])
    root=(theta-a)/(theta+a)
    check_s=[mp.mpf('.75')+1j*gamma,1+1j*gamma]
    return {
     "status":"all winding/Jensen/complex-source residuals NUMERICAL only",
     "height":mp.nstr(gamma,30),
     "scale":float(a),"contour_radius":float(r),
     "shift":float(theta),
     "fake_zero_local_w":mp.nstr(root,25),
     "fake_jensen_proved_lower_bound_given_critical_xi_zero":float(mp.log(r/abs(root))),
     "true_xi":sample(xi),"shifted_control":sample(fake),
     "complex_theta_source_abs_residuals":[mp.nstr(abs(theta_source(s)-xi(s)),8)
                                               for s in check_s],
     "positive_theta_kernel_phi_samples":[mp.nstr(phi_pos(x),9) for x in (0,mp.mpf('.5'),1)],
     "positive_shift_convolution_samples":[mp.nstr(conv_positive(u),9)
                                            for u in (0,mp.mpf('.5'))],
     "caveat":"No certified zero count; numerical winding is an experiment; complete local and global RH open."
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--numeric",action="store_true")
    parser.add_argument("--dps",type=int,default=38)
    parser.add_argument("--samples",type=int,default=160)
    args=parser.parse_args()
    if args.numeric:
        print(json.dumps(numerical_probe(args.dps,args.samples),indent=2))
    else:
        r0=abs(fake_root_coordinate(.25,.5))
        print(json.dumps({"fake_local_root":fake_root_coordinate(.25,.5),
                          "probe_radius":.45,
                          "Jensen_lower_bound":jensen_lower_bound(.45,r0),
                          "scope":"exact disk arithmetic; xi source test requires --numeric"},indent=2))
