"""Prime-side detector with explicit trivial-zero correction; finite zeros only a diagnostic."""
import contextlib,io,struct,json
from pathlib import Path
with contextlib.redirect_stdout(io.StringIO()):
 import detector_stationarity_audit as d
m=d.m
zeros=[]
for name in ['zeros_14.dat','zeros_5000.dat']:
 raw=Path(__file__).with_name(name).read_bytes();pos=8
 for _ in range(struct.unpack_from('<Q',raw)[0]):
  t0,t1,n0,n1=struct.unpack_from('<ddQQ',raw,pos);pos+=32;z=0
  assert n0==len(zeros)
  for j in range(n1-n0):
   z+=int.from_bytes(raw[pos:pos+13],'little');pos+=13
   zeros.append(m.mpf(t0)+m.mpf(z)*m.mpf(2)**-101)
 assert pos==len(raw)
M=3000
def R(u):
 return m.fsum(4*m.pi**4*(-m.expm1(-a))**2*m.exp(-a*(u-1))/(a*a*(a*a+4*m.pi*m.pi)**2) for a in [2*m.mpf(k)+m.mpf('.5') for k in range(1,M+1)])
# Decreasing a^-6: sum_(k>M)(2k+.5)^-6 <= integral_M^infinity (2x+.5)^-6 dx.
tail=4*m.pi**4/(10*(2*m.mpf(M)+m.mpf('.5'))**5)
rows=[]
for u in map(m.mpf,['1','1.5','2','3','5','7']):
 correction=R(u);prime=-d.K(u-1)-correction
 spectral=2*m.fsum(abs(d.S(1j*g))**2*m.cos(g*u) for g in zeros)
 rows.append(dict(u=str(u),R=str(correction),corrected_prime=str(prime),finite_zero_sum=str(spectral),difference=str(abs(prime-spectral))))
print(json.dumps(dict(scope=__doc__,zeros=len(zeros),zero_height=str(zeros[-1]),R_tail_upper_bound=str(tail),rows=rows),indent=2))
