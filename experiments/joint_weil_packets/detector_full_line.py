"""Full prime/gamma/pole kernel; zeros are used only for an independent diagnostic."""
import contextlib,io,json,struct
from pathlib import Path
with contextlib.redirect_stdout(io.StringIO()):
 import detector_stationarity_audit as d
m=d.m
m.mp.dps=50
G=d.G
c=d.S(m.mpf('.5'))*d.S(-m.mpf('.5'))
def C(u):
 u=abs(m.mpf(u));gu=G(u);L=u+2
 def integrand(t):
  if abs(t)<m.mpf('1e-22'): return -m.mpf('1.5')*gu
  return (2*m.exp(-2*t)*gu-m.exp(-t/2)*(G(u+t)+G(u-t)))/(-m.expm1(-2*t))
 knots=sorted(set([m.mpf(0),L]+[x for x in [u-1,1-u,u,u+1] if 0<x<L]))
 arch=-(m.log(m.pi)+m.euler)*gu+m.quad(integrand,knots)-gu*m.log1p(-m.exp(-2*L))
 limit=int(m.exp(u+1))
 assert limit<=d.N, 'Increase prime sieve; do not silently truncate.'
 prime=m.fsum(d.lam[n]/m.sqrt(n)*(G(m.log(n)-u)+G(m.log(n)+u)) for n in range(2,limit+1) if d.lam[n])
 return 2*c*m.cosh(u/2)-prime+arch
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
weights=[2*abs(d.S(1j*g))**2 for g in zeros]
rows=[]
cache={}
for u in map(m.mpf,['0','.2','.4','.6','.8','1','1.5','2']):
 v=C(u);cache[str(u)]=v
 z=m.fsum(w*m.cos(g*u) for w,g in zip(weights,zeros))
 rows.append(dict(u=str(u),prime_side=str(v),zero_diagnostic=str(z),difference=str(abs(v-z))))
A=m.matrix([[cache[str(m.mpf(abs(i-j))/5)] for j in range(5)] for i in range(5)])
eigs=m.eigsy(A,eigvals_only=True)
print(json.dumps(dict(scope=__doc__,rows=rows,toeplitz_eigenvalues=[str(x) for x in eigs],zeros=len(zeros)),indent=2))
