"""Pole-neutral three-box filter. Binary64 diagnostics, not an RH proof."""
import cmath
import json
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'odd_collision_trace'))
from audit import integral, prime_powers

a = math.log(2)
c = math.cosh(a/2)
weights = {0: 6.5, -1: -3*math.sqrt(2), 1: -3*math.sqrt(2), -2: 1, 2: 1}

def h(u):
    return max(0, 1-abs(u)/a)

def hphi(u):
    return math.fsum(w*h(u-k*a) for k, w in weights.items())

def H(z):
    return a if z == 0 else 2*(cmath.cosh(a*z)-1)/(a*z*z)

def Hphi(z):
    return 4*(cmath.cosh(a*z)-c)**2*H(z)


def main():
    checks=[]
    for z in (.5, -.5, .2+1.7j, .01+14.1j, 14.1j):
        direct=sum(integral(lambda u: (hphi(u)*cmath.exp(z*u)).real, k*a, (k+1)*a)
                   for k in range(-3, 3))
        direct+=1j*sum(integral(lambda u: (hphi(u)*cmath.exp(z*u)).imag, k*a, (k+1)*a)
                       for k in range(-3, 3))
        error=abs(direct-Hphi(z))
        assert error < 1e-9
        checks.append({'z':str(z), 'transform_error':error,'abs_transform':abs(Hphi(z))})
    pp=prime_powers(80000)
    def S(t):
        return math.fsum(w*h(y-t) for y,w in pp if abs(y-t)<=a)
    rows=[]
    for x in (10,100,1000,10000):
        t=math.log(x)
        five=math.fsum(w*S(t+k*a) for k,w in weights.items())
        direct=math.fsum(w*hphi(y-t) for y,w in pp if abs(y-t)<=3*a)
        assert abs(five-direct)<1e-9
        # D_phi(t)=integral h_phi(u-t) exp(-5u/2)/(1-exp(-2u)) du.
        # h_phi changes sign, but its exponential series has positive coefficients.
        tail=sum(integral(lambda v: hphi(v)*math.exp(-2.5*(t+v))/(-math.expm1(-2*(t+v))),
                          k*a,(k+1)*a) for k in range(-3,3))
        rows.append({'x':x,'five_scale_prime_sum':five,'D_phi':tail,
                     'C_phi':-five-tail,'identity_error':abs(five-direct)})
    z=.2+1.7j
    quartet=[z,z.conjugate(),-z,-z.conjugate()]
    corr=lambda t:sum(Hphi(w)*cmath.exp(w*t) for w in quartet).real
    maximum=max(corr(j/10) for j in range(601))
    assert maximum>100*max(1,abs(corr(0)))
    print(json.dumps({'status':'exploratory floating point; no all-scale estimate',
        'norm_squared_phi':6.5,'transform_checks':checks,'prime_checks':rows,
        'off_axis_control_C0':corr(0),'off_axis_control_max_to_60':maximum},indent=2))

if __name__=='__main__':
    main()
