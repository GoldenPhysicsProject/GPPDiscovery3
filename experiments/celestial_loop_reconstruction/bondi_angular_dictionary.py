"""Exact angular Bondi mass-constraint checks. Not a full spacetime solution."""
from fractions import Fraction as F
import json
from math import tanh, pi
from scipy.integrate import quad
# Coefficients in ascending powers; all angular identities checked over Q.
def add(a,b):
 c=[F(0)]*max(len(a),len(b))
 for i,v in enumerate(a):c[i]+=v
 for i,v in enumerate(b):c[i]+=v
 return c

def scale(a,s):return [s*v for v in a]
def mul(a,b):
 c=[F(0)]*(len(a)+len(b)-1)
 for i,v in enumerate(a):
  for j,w in enumerate(b):c[i+j]+=v*w
 return c

def der(a):return [i*a[i] for i in range(1,len(a))]
def equal(a,b):return all(v==0 for v in add(a,scale(b,-1)))
def integral(a):return sum(2*v/F(i+1) for i,v in enumerate(a) if i%2==0)
P0=[F(1)]; P2=[F(-1,2),0,F(3,2)]; P4=[F(3,8),0,F(-30,8),0,F(35,8)]
q=[F(1),0,F(-1)]
# T_AB=(D_A D_B - q_AB D^2/2)P2; T_theta theta=A.
A=scale(mul(q,der(der(P2))),F(1,2))
B=add(scale(mul(q,der(A)),-1),mul([0,F(2)],A))
divdiv=scale(der(B),-1)
norm=scale(mul(A,A),2)
assert equal(divdiv,scale(P2,12))
assert equal(norm,scale(mul(q,q),F(9,2)))
assert equal(mul(q,q),add(add(scale(P0,F(8,15)),scale(P2,F(-16,21))),scale(P4,F(8,35))))
assert integral(divdiv)==0 and integral(norm)==F(24,5)
# f=(1+tanh u)/2, F(u)=integral_{-infty}^u f'(v)^2 dv.
def f(u):return (1+tanh(u))/2
def fp(u):return (1-tanh(u)**2)/2
def flux(u):
 t=tanh(u);return (t-t**3/3+2/3)/4
rows=[]
for e in [.2,.7,1.]:
 # m=m0+3e f P2 -(9/16)e^2 flux (1-x^2)^2.
 # Delta m=-e^2/10 P0+(3e+e^2/7)P2-3e^2/70 P4.
 angular=quad(lambda x:4.5*(1-x*x)**2,-1,1,epsabs=1e-12)[0]*2*pi
 temporal=quad(lambda u:fp(u)**2,-20,20,epsabs=1e-12)[0]
 loss=-e*e*angular*temporal/(32*pi)
 assert abs(loss+e*e/10)<1e-12
 max_res=0
 for u in [-2,-.4,.3,1.5]:
  for x in [-.8,-.2,.4,.9]:
   p2=(3*x*x-1)/2
   def mass(v):return 2+3*e*f(v)*p2-9*e*e*flux(v)*(1-x*x)**2/16
   h=1e-4
   md=(-mass(u+2*h)+8*mass(u+h)-8*mass(u-h)+mass(u-2*h))/(12*h)
   rhs=3*e*fp(u)*p2-9*e*e*fp(u)**2*(1-x*x)**2/16
   max_res=max(max_res,abs(md-rhs))
 assert max_res<1e-9
 rows.append(dict(amplitude=e,mass_change=loss,expected=-e*e/10,max_finite_difference_constraint_residual=max_res))
results=dict(exact_angular_checks='passed over rational coefficients',sphere_norm_integral='48*pi/5',temporal_flux='1/3',tests=rows)
print(json.dumps(results,indent=2))
open('bondi_angular_dictionary_results.json','w').write(json.dumps(results,indent=2)+'\n')
