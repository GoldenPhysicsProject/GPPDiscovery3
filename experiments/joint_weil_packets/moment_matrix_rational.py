"""Exact Fraction audit of a conservative 17-mode tail-error bound.
L=4, T=1000, A<7, envelope b<4096, log(T)<10.
Analytic premises proved in companion note; not a Lean certificate.
"""
from fractions import Fraction as F
import json
T=1000;m=24;d=17;A=7
I=lambda p:F(4106,p-1)/T**(p-1)+F(1,(p-1)**2*T**(p-1))
rho=d*A**(2*m)*I(2*m+2)/(1-F(A,T))**2
trace=d*sum(A**(r+q)*I(r+q+2) for r in range(m) for q in range(m) if (r+q)%2==0)
delta=F(1,10**45)
checks={'delta_exceeds_rho':delta>rho,'squared_error_bound':(delta-rho)**2>4*rho*trace}
assert all(checks.values())
print(json.dumps({'scope':'Exact rational scalar inequality only; no finite-matrix PSD or infinite-complement certificate','m':m,'T':T,'dimension':d,'error_bound':'1e-45 before division by 2*pi','checks':checks,'rho':{'numerator':str(rho.numerator),'denominator':str(rho.denominator)},'trace_upper':{'numerator':str(trace.numerator),'denominator':str(trace.denominator)}},indent=2))
