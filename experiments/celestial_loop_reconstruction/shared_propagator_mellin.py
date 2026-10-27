"""Damped fixed-angle null-line propagator kernel, NOT a full loop amplitude.
M denotes mass squared. Mellin powers are energy^(alpha-1).
"""
import json
import mpmath as mp
import numpy as np
from scipy.special import roots_laguerre
mp.mp.dps=25

def direct(alpha,beta,a,b,M,c):
    # Analytically integrate y using Tricomi U; outer energy quadrature.
    def f(x):
        inner=mp.gamma(beta)*(M/(c*x))**beta/M*mp.hyperu(beta,beta,b*M/(c*x))
        return x**(alpha-1)*mp.exp(-a*x)*inner
    return mp.quad(f,[0,1,mp.inf])

def spectral(alpha,beta,a,b,M,c):
    sigma=mp.mpf('0.25')
    def f(t):
        s=sigma+1j*t
        return (mp.gamma(s)*mp.gamma(1-s)*M**(s-1)*c**(-s)
                *mp.gamma(alpha-s)*a**(s-alpha)
                *mp.gamma(beta-s)*b**(s-beta))/(2*mp.pi)
    return mp.quad(f,[-mp.inf,-10,-2,0,2,10,mp.inf])

rows=[]
for lam,nu in [(0,0),(.7,-.4),(1.2,-1.2)]:
    alpha=1+1j*mp.mpf(str(lam)); beta=1+1j*mp.mpf(str(nu))
    args=(alpha,beta,mp.mpf('1.3'),mp.mpf('.8'),mp.mpf('1.1'),mp.mpf('.6'))
    d=direct(*args); s=spectral(*args)
    rel=abs(d-s)/abs(d)
    assert rel<mp.mpf('1e-18')
    rows.append({'lambda':lam,'nu':nu,'energy':str(d),'mellin_barnes':str(s),'relative_error':str(rel)})

# Independent 2D energy quadrature, no hyperu or Mellin-Barnes.
pair=[]
chain=[]
for n in (32,64,128,256):
    x,w=roots_laguerre(n)
    pair.append({'nodes':n,'value':float(np.einsum('i,j,ij->',w,w,1/(1.1+.6*x[:,None]*x[None,:]/(1.3*.8)))/(1.3*.8))})
    # Chain of N energies, alpha_i=1, a_i=1, all M_e=1.1,c_e=.6.
    K=1/(1.1+.6*x[:,None]*x[None,:])
    f=np.ones(n); vals={1:float(w@f)}
    for N in range(2,11):
        f=K@(w*f)
        vals[N]=float(w@f)
    chain.append({'nodes':n,'chain_values':vals})
print(json.dumps({'scope':'fixed-angle damped cut-energy kernel; chain length is NOT loop count',
                 'pair_tests':rows,'independent_pair_quadrature':pair,'chain_tests':chain},indent=2))
