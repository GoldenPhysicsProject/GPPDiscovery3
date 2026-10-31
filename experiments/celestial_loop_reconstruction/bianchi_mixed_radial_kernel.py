"""Numerical cross-check only; no cusp coefficients or certified intervals."""
import json
import mpmath as mp
mp.mp.dps = 40

def kernel(s, r, t, a):
    return (2**(s-3)*a**(-s)/mp.gamma(s) *
            mp.fprod(mp.gamma((s+e*1j*r+f*1j*t)/2)
                     for e in (-1,1) for f in (-1,1)))

def check(s, r, t, a):
    # Independent Bessel quadrature, rescaled x=a*y.
    q = a**(-s)*mp.quad(lambda x: x**(s-1)*mp.besselk(1j*r,x)*mp.besselk(1j*t,x),
                       [0, mp.mpf('.1'), 1, 5, 15, 40, mp.inf])
    g = kernel(s,r,t,a)
    err=abs(q-g)/abs(g)
    assert err < mp.mpf('1e-22'), err
    return {k:str(v) for k,v in dict(s=s,r=r,t=t,a=a,quadrature=q,gamma_formula=g,relative_error=err).items()}

if __name__ == '__main__':
    r0=mp.mpf('6.6221193402528') # Then's approximate numerical eigenparameter, not certified.
    rows=[check(mp.mpf(3),mp.mpf('.7'),mp.mpf('1.2'),2*mp.pi),
          check(mp.mpf(3),r0,r0,2*mp.pi),
          check(2+mp.mpf('.4')*1j,mp.mpf('.7'),mp.mpf('1.2'),mp.mpf('1.3'))]
    print(json.dumps({'dps':mp.mp.dps,'scope':'radial kernel only; not a full automorphic overlap', 'checks':rows},indent=2))
