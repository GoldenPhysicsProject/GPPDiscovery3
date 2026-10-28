import numpy as np
from scipy.integrate import quad
import json
pi=np.pi
b=.73
def F(t): return pi**(-.25)*np.exp(-(t-b)**2/2)
rows=[]
for d in (1,2,3):
    norm=quad(lambda l: abs(F(l/d)/np.sqrt(d))**2,-np.inf,np.inf)[0]
    errors=[]
    for l in (-4.,-.3,0.,2.7):
        t=l/d
        for p in (2,3,5,11):
            a=np.log(p)
            # Independent integration of translated log-coordinate Gaussian.
            def integrand(u):
                v=u-a
                return pi**(-.25)*np.exp(-v*v/2+1j*b*v-1j*t*u)/np.sqrt(2*pi*d)
            value=quad(lambda u: integrand(u).real,-15,15,epsabs=1e-11)[0]+1j*quad(lambda u: integrand(u).imag,-15,15,epsabs=1e-11)[0]
            expected=np.exp(-1j*t*a)*F(t)/np.sqrt(d)
            errors.append(abs(value-expected))
    rows.append({'d':d,'norm':norm,'max_prime_translation_error':max(errors)})
a=1.
overlap=quad(lambda x: np.exp(-(x*x+(x-a)**2)/2)/np.sqrt(pi),-np.inf,np.inf)[0]
result={'rows':rows,'angular_translation_overlap':overlap,'angular_range_leakage_squared':1-overlap**2,'expected_leakage_squared':1-np.exp(-a*a/2)}
assert all(abs(r['norm']-1)<1e-10 and r['max_prime_translation_error']<1e-10 for r in rows)
assert abs(overlap-np.exp(-a*a/4))<1e-12
print(json.dumps(result,indent=2))
