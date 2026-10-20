"""Signed decomposition and deterministic centering of prime Cesaro energy. Diagnostic only."""
import math,json
import mpmath as mp
mp.mp.dps=70
N=100000
sieve=bytearray(b'\1')*(N+1);sieve[:2]=b'\0\0'
for p in range(2,math.isqrt(N)+1):
 if sieve[p]:sieve[p*p:N+1:p]=b'\0'*((N-p*p)//p+1)
lam={}
for p in range(2,N+1):
 if sieve[p]:
  n=p
  while n<=N:lam[n]=mp.log(p);n*=p
cuts={10,100,1000,10000,100000}
points=sorted({1,N}|set(lam)|cuts)
psi=mp.mpf(0);weighted=mp.mpf(0);energy=mp.mpf(0);checks=[];rows=[]
PP=mp.mpf(0);PC=mp.mpf(0);CC=mp.mpf(0);meanJ=mp.mpf(0)
kappa=2+mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi)
for left,right in zip(points,points[1:]):
 w=lam.get(left,mp.mpf(0));psi+=w;weighted+=w/mp.sqrt(left)
 A=2*weighted+2;B=-2*psi;C=mp.mpf(-2)
 a=mp.log(left);d=mp.log(mp.mpf(right)/left)
 def integ(k):
  return d if k==0 else mp.exp(k*a)*mp.expm1(k*d)/k
 piece=(A*A+2*B*C)*d+B*B*integ(-1)+C*C*integ(1)+2*A*B*integ(-mp.mpf('.5'))+2*A*C*integ(mp.mpf('.5'))
 assert piece>=0
 if right<=100:
  direct=mp.quad(lambda t:(A+B*mp.exp(-t/2)+C*mp.exp(t/2))**2,[a,a+d])
  checks.append(abs(piece-direct));assert abs(piece-direct)<mp.mpf('1e-60')
 energy+=piece
 P0=A-2
 PP+=P0*P0*d+B*B*integ(-1)+2*P0*B*integ(-mp.mpf('.5'))
 PC+=2*P0*integ(mp.mpf('.5'))-2*P0*d+2*B*d-2*B*integ(-mp.mpf('.5'))
 CC+=4*integ(1)-8*integ(mp.mpf('.5'))+4*d
 meanJ+=A*d+B*integ(-mp.mpf('.5'))+C*integ(mp.mpf('.5'))
 if right in cuts:
  t=mp.log(right);J=A+B/mp.sqrt(right)+C*mp.sqrt(right)
  assert abs(energy-(PP-2*PC+CC))<mp.mpf('1e-55')
  centered=energy-2*kappa*meanJ+kappa*kappa*t
  assert centered>=0
  rows.append(dict(prime_prime=mp.nstr(PP,30),twice_prime_continuum=mp.nstr(2*PC,30),continuum_continuum=mp.nstr(CC,30),centered_energy=mp.nstr(centered,30),raw_to_unsigned=mp.nstr(energy/(PP+2*PC+CC),15),X=right,T=mp.nstr(t,20),J=mp.nstr(J,30),integrated_J_squared=mp.nstr(energy,30),cesaro_mean_square=mp.nstr(energy/t,30)))
print(json.dumps(dict(scope='Finite prime-only numerical diagnostics; not certified and not a uniform bound',kappa=mp.nstr(kappa,40),precision=mp.mp.dps,prime_power_count=len(lam),quadrature_checks=len(checks),max_check_error=mp.nstr(max(checks),6),rows=rows),indent=2))
