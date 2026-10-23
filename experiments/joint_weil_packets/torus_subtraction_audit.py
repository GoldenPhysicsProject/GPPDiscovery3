"""Dual Poisson evaluation with analytic truncation bound (not interval arithmetic)."""
import mpmath as m,itertools,json
m.mp.dps=60
rows=[]
for N in [2,3]:
 w=[m.log(p) for p in [2,3,5][:N]];r2=sum(a*a for a in w);L=10;t=m.mpf(1)
 lattice=[]
 for v in itertools.product(range(-L,L+1),repeat=N):
  if not any(v):continue
  norm=sum(k*k for k in v);dot=sum(k*a for k,a in zip(v,w));parallel=dot*dot/r2
  lattice.append((m.mpf(norm)-parallel,parallel))
 for e in map(m.mpf,['.1','.01','.001','.0001']):
  V=(m.pi/t)**(m.mpf(N)/2)/e**(m.mpf(N-1)/2)/m.sqrt(r2)
  B=V/m.sqrt(1+e/r2)
  rem=B*m.fsum(m.exp(-m.pi**2/t*(perp/e+par/(e+r2))) for perp,par in lattice)
  c=m.pi**2/(t*(e+r2))
  tail=B*N*(1+m.sqrt(m.pi/c))**(N-1)*m.sqrt(m.pi/c)*m.erfc(m.sqrt(c)*L)
  # Stable exact zero-mode minus leading asymptotic subtraction.
  bgdiff=V*m.expm1(-m.log1p(e/r2)/2)
  rows.append(dict(N=N,epsilon=str(e),Z_minus_leading=str(bgdiff+rem),Z_minus_exact_background=str(rem),dual_tail_bound=str(tail)))
print(json.dumps(dict(scope=__doc__,rows=rows),indent=2))
