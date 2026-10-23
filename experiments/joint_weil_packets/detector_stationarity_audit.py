import mpmath as m, math, json
m.mp.dps=40
def S(z):
 if z==0:return m.mpf('.5')
 return 2*m.pi**2*m.expm1(z)/(z*(z*z+4*m.pi**2))
def G(v):
 v=abs(v)
 if v>=1:return m.mpf(0)
 return (1-v)/4+(1-v)*m.cos(2*m.pi*v)/8+3*m.sin(2*m.pi*v)/(16*m.pi)
for v in ['0','.2','.7']:
 v=m.mpf(v)
 direct=m.quad(lambda x:m.sin(m.pi*x)**2*m.sin(m.pi*(x+v))**2,[0,1-v])
 assert abs(G(v)-direct)<m.mpf('1e-35')
N=math.ceil(math.exp(8))
lam=[m.mpf(0)]*(N+1)
for p in range(2,N+1):
 if all(p%d for d in range(2,math.isqrt(p)+1)):
  k=p
  while k<=N:
   lam[k]=m.log(p);k*=p
c=S(m.mpf('.5'))*S(m.mpf('-.5'))
def K(t):
 return sum(lam[n]/m.sqrt(n)*G(m.log(n)-t-1) for n in range(2, min(N,int(m.exp(t+2)))+1) if lam[n])-m.exp((t+1)/2)*c
rows=[dict(t=str(t),K=str(K(m.mpf(t)))) for t in ['0','.1','.25','.5','1','2','3','4','5','6']]
print(json.dumps(dict(rows=rows,K0=str(K(0)),scope='Exact finite prime sum evaluated numerically; no RH assumed'),indent=2))
