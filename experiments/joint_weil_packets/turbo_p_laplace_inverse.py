"""High-precision finite dictionary and inverse checks; no RH bound."""
import mpmath as m,json
m.mp.dps=60
c=3-1/m.expm1(m.mpf('.5'));b=m.mpf('.5')+1/c
rows=[]
for z in [m.mpc('0.9','.3'),m.mpc('1.7','2'),m.mpc('2.1','0')]:
 P=(1-c*(z-m.mpf('.5')))/((z-m.mpf('.5'))**2*m.zeta(z+m.mpf('.5')))-1/(z-m.mpf('.5'))
 inverse=-P/c+P/((c+1)*(z+m.mpf('.5')))+(-c-P)/(c*c*(c+1)*(z-b))-1/((c+1)*(z+m.mpf('.5')))
 exact=1/((z+m.mpf('.5'))*m.zeta(z+m.mpf('.5')))
 assert abs(inverse-exact)<m.mpf('1e-50')
 rows.append({'z':str(z),'inverse_transform_error':str(abs(inverse-exact))})
N=10000;mu=[1]*(N+1);mu[0]=0;prime=[True]*(N+1);prime[0]=prime[1]=False
for p in range(2,N+1):
 if prime[p]:
  for n in range(p,N+1,p):
   mu[n]*=-1
   if n>p:prime[n]=False
  for n in range(p*p,N+1,p*p):mu[n]=0
# p(u)=exp(u/2)*(B_n*u-c*B_n-L_n-1) on log n <=u<log(n+1).
def primitive(u,k,B,C):
 return m.exp(k*u)*(B*(u/k-1/k**2)+C/k)
B=m.mpf(0);L=m.mpf(0);M=0;past=m.mpf(0);weighted=m.mpf(0);samples=[]
for n in range(1,N+1):
 B+=m.mpf(mu[n])/n;L+=m.mpf(mu[n])*m.log(n)/n;M+=mu[n]
 y=m.log(n);P=m.sqrt(n)*((y-c)*B-L-1)
 if n in [1,10,100,1000,10000]:
  causal=m.exp(-y/2)*past
  future=m.exp(b*y)*(-c-weighted) # exact total moment LP(b)=-c
  recovered=-P/c+causal/(c+1)+future/(c*c*(c+1))-m.exp(-y/2)/(c+1)
  actual=M/m.sqrt(n)
  assert abs(recovered-actual)<m.mpf('1e-45')
  samples.append({'X':n,'p':str(P),'M_normalized':str(actual),'inverse_error':str(abs(recovered-actual))})
 if n<N:
  u=m.log(n+1);C=-c*B-L-1
  past+=primitive(u,m.mpf(1),B,C)-primitive(y,m.mpf(1),B,C)
  k=m.mpf('.5')-b
  weighted+=primitive(u,k,B,C)-primitive(y,k,B,C)
print(json.dumps({'scope':__doc__,'c':str(c),'b':str(b),'transform_checks':rows,'arithmetic_checks':samples},indent=2))
