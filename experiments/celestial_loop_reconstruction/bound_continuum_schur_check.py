"""Friedrichs-model control, not automorphic coefficient data or interval certification."""
import json
import mpmath as m
m.mp.dps=40
C=lambda x:x*m.exp(-x)
A=lambda x:1+x*x
S=lambda z:m.quad(lambda x:abs(C(x))**2/(A(x)-z),[0,1,4,m.inf])
s0=S(0); energy=m.mpf(2); gc=m.sqrt(energy/s0)
rows=[]
for ratio in [m.mpf('.5'),m.mpf('1'),m.mpf('2')]:
    g=ratio*gc
    schur=energy-g*g*s0
    h=lambda x:-g*C(x)/A(x)
    direct=energy+m.quad(lambda x:A(x)*abs(h(x))**2+2*g*C(x)*h(x),[0,1,4,m.inf])
    assert abs(direct-schur)<m.mpf('1e-35')
    rows.append({'g_over_gcrit':str(ratio),'schur':str(schur),'square_completion_error':str(abs(direct-schur))})
g=2*gc
root=m.findroot(lambda z:energy-z-g*g*S(z),(-8,-1))
assert root<0
# Same modulus and self energy under a nontrivial symmetric unitary phase.
B=lambda x,a,gamma:((x-gamma-1j*a)/(x-gamma+1j*a))*((x+gamma-1j*a)/(x+gamma+1j*a))
phase_errors=[]
for a in [m.mpf('.5'),m.mpf('.7')]:
    sp=m.quad(lambda x:abs(B(x,a,m.mpf(3))*C(x))**2/A(x),[0,1,4,m.inf])
    phase_errors.append(str(abs(sp-s0)))
    assert abs(sp-s0)<m.mpf('1e-35')
print(json.dumps({'scope':'specified toy form factor C=lambda*exp(-lambda), not Gaussian cusp data','dps':m.mp.dps,'sigma0':str(s0),'gcrit':str(gc),'checks':rows,'negative_eigenvalue_at_twice_gcrit':str(root),'phase_invariance_errors':phase_errors},indent=2))
