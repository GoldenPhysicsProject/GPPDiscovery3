"""Independent finite checks of the balanced prime-error integral.
Numerical diagnostics only; analytic bounds are proved in the companion note.
"""
import json, math
import mpmath as mp
mp.mp.dps=45
N=1000
prime=[True]*(N+1);prime[0]=prime[1]=False
lam=[mp.mpf(0)]*(N+1)
for p in range(2,N+1):
 if prime[p]:
  for n in range(p*p,N+1,p):prime[n]=False
  n=p
  while n<=N:lam[n]=mp.log(p);n*=p

def power_integral(a,b,p):
 return mp.log(b/a) if p==0 else (b**p-a**p)/p

def error_integral(X,p):
 # integral_1^X (psi(x)-x) x^(p-1) dx, integrated cell by cell
 total=mp.mpf(0);psi=mp.mpf(0)
 for n in range(1,int(mp.ceil(X))):
  psi+=lam[n];a=mp.mpf(n);b=min(a+1,X)
  total+=psi*power_integral(a,b,p)-power_integral(a,b,p+1)
 return total
rows=[]
for X in map(mp.mpf,['10','100','1000']):
 T=mp.log(X)
 for eps in [mp.mpf(0),mp.mpf('.25'),mp.mpf('.5'),mp.mpf('.5')-mp.mpf('.4')/mp.sqrt(T)]:
  direct=2*mp.fsum(lam[n]*(mp.mpf(n)**(-mp.mpf('.5')-eps)-X**(-mp.mpf('.5'))*mp.mpf(n)**(-eps)) for n in range(1,int(X)+1))
  counter=2*(power_integral(mp.mpf(1),X,mp.mpf('.5')-eps)-X**(-mp.mpf('.5'))*power_integral(mp.mpf(1),X,1-eps))+2*(1-X**(-mp.mpf('.5')))
  direct-=counter
  integral=(1+2*eps)*error_integral(X,-mp.mpf('.5')-eps)-2*eps*X**(-mp.mpf('.5'))*error_integral(X,-eps)
  error=abs(direct-integral);assert error<mp.mpf('1e-37')
  rows.append({'X':str(X),'epsilon':str(eps),'J':str(direct),'identity_error':str(error)})
A=mp.mpf('9.39');b=mp.mpf('1.515');c=mp.mpf('.8274')
# A tighter non-certified numerical value for the analytic uniform constant.
B_tight=2*(1+2*A*mp.gamma(2*b+2)/(c-mp.mpf('.4'))**(2*b+2))
# The rational constant in the note uses A<10, u^b<=1+u^2, c-.4>.4.
from fractions import Fraction as Q
B=2*(1+10*(2/Q(2,5)**2+240/Q(2,5)**6))
assert B==1172127
print(json.dumps({'scope':__doc__,'rows':rows,'exact_uniform_bound':str(B),'tighter_constant_numeric_only':str(B_tight),'endpoint_limit':str(-2*(1+mp.euler))},indent=2))
