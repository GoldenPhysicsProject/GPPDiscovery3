# Astra recovery audit, 2026-10-27.
# Source: Daniel's attached 15-verify_loops_from_cuts.py, preserved below.
# Caveats: section 5 compares definitions with themselves, not independent pole tests.
# ROOT/subprocess are unused: this script does NOT execute the original step scripts.
# Numerical checks are not interval certificates or Lean proofs.
"""Independent high-precision verification suite for the celestial cut/dispersion paper.

This is a numerical/symbolic verification suite, not a formal proof.  It checks
exact identities, the celestial phase-space normalization, the Mellin beta
identity, the regulated cut, the Mellin pole coefficients, and the dispersion
reconstruction.  The original step scripts are also executed as cross-checks.
"""
import subprocess, sys
from pathlib import Path
import sympy as sp
import mpmath as mp
mp.mp.dps = 35
ROOT = Path(__file__).resolve().parent.parent / "unpack18"

def close(a,b,tol=mp.mpf('1e-35')):
    return abs(a-b) <= tol*max(1,abs(a),abs(b))

# 1. Null celestial parametrization and one-particle Jacobian.
x,y,w=sp.symbols('x y w', positive=True)
q=sp.Matrix([1+x*x+y*y,2*x,2*y,1-x*x-y*y])
assert sp.simplify(q.dot(sp.diag(1,-1,-1,-1)*q)) == 0
J=sp.Matrix([[sp.diff(w*q[i],v) for v in (w,x,y)] for i in (1,2,3)]).det()
assert sp.simplify(J/(2*w*q[0]) - 2*w) == 0

# 2. Round-sphere phase-space normalization.
I = mp.quad(lambda r: (2*mp.pi*r)/(8*mp.pi**2*(1+r*r)**2), [0,mp.inf])
assert close(I,1/(8*mp.pi))

# 3. Mellin phase-space identity, independently via Euler beta integral.
D5=mp.mpc('1.27','0.41'); D6=mp.mpc('0.83','-0.19'); M=mp.mpf('3.7')
closed=(1/(8*mp.pi))*(M/2)**(D5+D6-2)*mp.gamma(D5)*mp.gamma(D6)/mp.gamma(D5+D6)
beta=(1/(8*mp.pi))*(M/2)**(D5+D6-2)*mp.quad(lambda u:u**(D5-1)*(1-u)**(D6-1),[0,1])
assert close(closed,beta,mp.mpf('1e-25'))

# 4. Regulated cut: celestial/Feynman-parameter closed form agrees with direct x integral.
kappa=mp.mpf('1.73'); mu2=mp.mpf('0.037'); s=mp.mpf('2.41')
c=mu2*(s+mu2); d=s*s*kappa
cut_int=(1/(8*mp.pi))*mp.quad(lambda u:1/(c+d*u*(1-u)),[0,1])
r=mp.sqrt(d/(d+4*c))
cut_closed=(1/(8*mp.pi))*(4/mp.sqrt(d*(d+4*c)))*mp.atanh(r)
assert close(cut_int,cut_closed,mp.mpf('1e-25'))

# 5. Universal Mellin pole coefficients from the exact asymptotic tail.
L=mp.log(kappa/mu2)
A2=1/(4*mp.pi*kappa)
A1=L/(4*mp.pi*kappa)
assert close(A2,1/(4*mp.pi*kappa))
assert close(A1,L/(4*mp.pi*kappa))

# 6. Cutkosky normalization anchored on the massless bubble.
# Bubble Im J_2 = pi and C_2 = 1/(8 pi), hence Im J_2/C_2 = 8 pi^2.
assert close(mp.pi/(1/(8*mp.pi)),8*mp.pi**2)

# 7. Independent dispersive reconstruction against the direct Feynman-parameter box.
def G(c,d):
    w0=d+4*c
    return 8*mp.atanh(mp.sqrt(d/w0))/(mp.sqrt(d)*w0**mp.mpf('1.5')) + 2/(c*w0)

def J_direct(s,u,mu2):
    S,U=-s,-u
    def inner(x2,x4):
        rho=1-x2-x4
        if rho<=0: return mp.mpf('0')
        c0=mu2*(x2+x4)+U*x2*x4
        d0=S*rho**2
        return rho*G(c0,d0)
    return mp.quad(lambda x2: mp.quad(lambda v:(1-x2)*inner(x2,(1-x2)*v),[0,.5,1]),[0,.01,.5,.99,1])

def C_cut(sp,u,mu2):
    c0=mu2*(sp+mu2); d0=sp*(-u); w0=d0+4*c0
    return (1/(8*mp.pi))*(4/mp.sqrt(d0*w0))*mp.atanh(mp.sqrt(d0/w0))

def J_disp(s,u,mu2):
    S=-s
    return 8*mp.pi*mp.quad(lambda sp:C_cut(sp,u,mu2)/(sp+S),[0,mu2,1,10,100,1000,mp.inf])

for s0,u0,m20 in [(-3,-2,.5)]:
    s0,u0,m20=map(mp.mpf,(s0,u0,m20))
    jd=J_direct(s0,u0,m20); jp=J_disp(s0,u0,m20)
    assert close(jd,jp,mp.mpf('1e-12'))

print('PASS: independent Feynman-parameter box = Cutkosky dispersion reconstruction at the tested Euclidean point')
print('PASS: celestial Jacobian and null geometry')
print('PASS: two-particle celestial phase-space normalization = 1/(8*pi)')
print('PASS: Mellin phase-space Gamma/Beta identity')
print('PASS: regulated box cut closed form')
print('PASS: universal Mellin double/simple-pole coefficients')
print('PASS: Cutkosky constant 8*pi^2 anchored on the bubble')
