"""Exactly designed Friedrichs control; not an arithmetic cusp coefficient."""
import json
import mpmath as m
m.mp.dps=40
g=m.mpf(1); E=m.mpf(2); mu=E-g*g/4
A=lambda x:1+x*x
C=lambda x:(x*x-1)*m.exp(-x)
Z=1/(1+g*g/2); a=m.sqrt(Z)
f=lambda x:-g*a*m.exp(-x)
I=m.quad(lambda x:(x*x-1)*m.exp(-2*x),[0,1,m.inf])
norm=a*a+m.quad(lambda x:abs(f(x))**2,[0,1,m.inf])
bound_res=mu*a+g*m.quad(lambda x:C(x)*f(x),[0,1,m.inf])-E*a
cont_res=max(abs(A(x)*f(x)+g*a*C(x)-E*f(x)) for x in [m.mpf('.1'),m.mpf(1),m.mpf(2),m.mpf(10)])
s0=m.quad(lambda x:C(x)**2/A(x),[0,1,m.inf])
assert abs(I+m.mpf(1)/4)<m.mpf('1e-35')
assert abs(norm-1)<m.mpf('1e-35')
assert max(abs(bound_res),cont_res)<m.mpf('1e-35')
assert mu-g*g*s0>0
print(json.dumps({k:str(v) for k,v in dict(scope='engineered control, not a Gaussian cusp mode',mu=mu,E=E,Z=Z,norm=norm,Sigma_E=I,bound_residual=abs(bound_res),continuum_residual=cont_res,positivity_scalar=mu-g*g*s0,survival_probability_lower_bound=(2*Z-1)**2).items()},indent=2))
