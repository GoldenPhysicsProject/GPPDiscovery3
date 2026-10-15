#!/usr/bin/env python3
"""Exact rational certificate for a real-test chi7 cubic fixed-window bound.
Astra / GPPDiscovery3. No floating-point branch decisions; no RH/GRH claim.
Analytic ground-state and block reduction are documented in companion note.
"""
from fractions import Fraction as F
class I:
 def __init__(self,l,u=None):
  lo=F(l);hi=F(l if u is None else u);scale=10**25
  self.l=F((lo.numerator*scale)//lo.denominator,scale)
  self.u=F(-((-hi.numerator*scale)//hi.denominator),scale)
 def __add__(self,o):
  o=asI(o);return I(self.l+o.l,self.u+o.u)
 __radd__=__add__
 def __neg__(self): return I(-self.u,-self.l)
 def __sub__(self,o):return self+-asI(o)
 def __rsub__(self,o):return asI(o)+-self
 def __mul__(self,o):
  o=asI(o);v=[a*b for a in (self.l,self.u) for b in (o.l,o.u)];return I(min(v),max(v))
 __rmul__=__mul__
 def __truediv__(self,o):
  o=asI(o);assert o.l>0 or o.u<0;return self*I(1/o.u,1/o.l)
 def __pow__(self,n):
  assert n>=0;r=I(1)
  for _ in range(n):r=r*self
  return r
 def show(self):return (float(self.l),float(self.u))
def asI(x):return x if isinstance(x,I) else I(x)
def ath(x,n=180):
 x=F(x);assert 0<=x<1
 s=sum((x**(2*j+1)/ (2*j+1) for j in range(n)),F(0))
 return I(s,s+x**(2*n+1)/((2*n+1)*(1-x*x)))
def at(x,n=180):
 x=F(x);assert 0<=x<1
 s=sum(((-1)**j*x**(2*j+1)/(2*j+1) for j in range(n)),F(0))
 nxt=(-1)**n*x**(2*n+1)/(2*n+1)
 return I(min(s,s+nxt),max(s,s+nxt))
def logq(x):
 x=F(x);assert x>=1
 # reduce to [1,2] so convergence remains fast
 m=0
 while x>2:x/=2;m+=1
 return 2*ath((x-1)/(x+1))+2*m*ath(F(1,3))
def mono(fn,x):return I(fn(x.l).l,fn(x.u).u)
pi=16*at(F(1,5))-4*at(F(1,239))
h=logq(3)/2;t=logq(2)
sqrt2=I('1.41421356237309504880','1.41421356237309504881')
sqrt3=I('1.73205080756887729352','1.73205080756887729353')
assert sqrt2.l**2<2<sqrt2.u**2 and sqrt3.l**2<3<sqrt3.u**2
z=I('0.75983568565159254733','0.75983568565159254734')
assert 3*z.l**4<1<3*z.u**4
T=2*(mono(ath,z)+mono(at,z))
# gamma <= H_n-log n, strict upper bound; lower endpoint not needed.
gamma_upper=(I(sum((F(1,j) for j in range(1,1001)),F(0)))-logq(1000)).u
c_upper=(mono(logq,8*pi)+gamma_upper+pi/2).u
d=T+logq(7)-c_upper
alpha=9*h/(4*sqrt3) # 2h*J_even(2h), since exp(-2h)=1/3
b=t/(2*sqrt2) # coefficient norm of B=(kappa/2)(T_shift+T_shift*)
r=(2*h-t)/h
beta=b*r
v=b*b*r-beta*beta # squared off-diagonal block norm
q=F(1,200)
A=d+beta-q;D=d+alpha-b-q
assert A.l>0 and D.l>0 and (A*D-v).l>0
for name,obj in [('h',h),('d (gamma upper used)',d),('alpha',alpha),('b',b),('beta',beta),('offdiag_squared',v),('A',A),('D',D),('determinant',A*D-v)]:print(name,obj.show())
print('PASS: exact rational comparisons give Q_chi7(g) > (1/200)||g||^2.')
print('Real smooth tests, both primitive cubic mod-7 characters, h <= log(3)/2.')
print('Analytic reduction requires independent audit; not Lean-verified; RH/GRH OPEN.')
