"""Free radiative Mellin reconstruction tests; not interacting quantum gravity."""
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma

def cq(f,a,b):
    return quad(lambda x:f(x).real,a,b,epsabs=2e-10,limit=400)[0]+1j*quad(lambda x:f(x).imag,a,b,epsabs=2e-10,limit=400)[0]
rows=[]
for a,u in [(1.,.3),(.7,-1.2),(2.,0.)]:
    # f(omega)=exp(-a omega); F(Delta)=Gamma(Delta)/a**Delta.
    # epsilon is an explicit plane-wave convergence regulator.
    eps=.2
    def integrand(lam):
        d=1+1j*lam
        return gamma(d)*a**(-d)*gamma(2-d)*(eps+1j*u)**(d-2)/(2*np.pi)
    spectral=cq(integrand,-70,70)
    momentum=cq(lambda w:w*np.exp(-(a+eps+1j*u)*w),0,np.inf)
    exact=1/(a+eps+1j*u)**2
    assert abs(spectral-exact)<1e-9 and abs(momentum-exact)<1e-9
    rows.append(dict(a=a,u=u,epsilon=eps,spectral=[spectral.real,spectral.imag],momentum=[momentum.real,momentum.imag],error=abs(spectral-exact)))
print(json.dumps(rows,indent=2))
open('celestial_dictionary_check_results.json','w').write(json.dumps(rows,indent=2)+'\n')
