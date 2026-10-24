"""Theta half-wave identities and endpoint cancellation; floating diagnostics only."""
import json
import mpmath as m
m.mp.dps=60

def psi(x):return m.fsum(m.exp(-m.pi*n*n*x) for n in range(1,9))
def phi(u):
    return m.fsum((2*m.pi**2*n**4*m.exp(9*u)-3*m.pi*n*n*m.exp(5*u))*m.exp(-m.pi*n*n*m.exp(4*u)) for n in range(1,9))
def A(z):return 4*m.quad(lambda u:phi(u)*m.exp(2*z*u),[0,.25,.5,.75,1,2])
def I(s):return 4*m.quad(lambda u:psi(m.exp(4*u))*m.exp(2*s*u),[0,.25,.5,.75,1,2])
def xi(s):return s*(s-1)*m.pi**(-s/2)*m.gamma(s/2)*m.zeta(s)/2
rows=[]
for s in [m.mpc('.5','5'),m.mpc('.6','14'),m.mpc('.3','21')]:
    z=s-m.mpf('.5');a=A(z)
    rhs=m.mpf('.25')+z*psi(1)+s*(s-1)*I(s)/2
    err=abs(a-rhs);sumerr=abs(a+A(-z)-xi(s))
    assert err<m.mpf('1e-50') and sumerr<m.mpf('1e-50')
    rows.append(dict(s=str(s),half_identity_error=str(err),xi_identity_error=str(sumerr)))
phi0=phi(0);dd=m.diff(phi,0,2);asym=[]
for t in [20,40,80]:
    a=A(1j*t);x=xi(m.mpf('.5')+1j*t)
    leading=2j*phi0/t
    nextterm=leading-1j*dd/(2*t**3)
    asym.append(dict(t=t,A_imag=str(a.imag),Xi_real=str(x.real),leading_relative_error=str(abs((a-leading)/a)),two_term_relative_error=str(abs((a-nextterm)/a)),cancellation_ratio=str(2*abs(a)/abs(x))))
print(json.dumps(dict(scope=__doc__,phi_zero=str(phi0),identity_checks=rows,asymptotic_checks=asym),indent=2))
