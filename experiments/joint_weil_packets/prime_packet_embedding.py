"""Diagnostic quadratures for the exact packet construction, not an RH test."""
import mpmath as m, json
m.mp.dps=60
rows=[]
for eps in [m.mpf('1'),m.mpf('.01')]:
 for n in [1,7,100]:
  c=m.log(n); h=eps*m.log1p(m.mpf(1)/n)/4
  mean=m.quad(lambda y:c+y,[-h,h])/(2*h)
  var=m.quad(lambda y:y*y,[-h,h])/(2*h)
  assert abs(mean-c)<m.mpf('1e-55')
  assert abs(var-h*h/3)<m.mpf('1e-55')
  t=m.mpf('2.3')
  amp=m.quad(lambda y:m.exp(1j*t*(c+y)),[-h,h])/(2*h)
  exact=m.exp(1j*t*c)*m.sin(t*h)/(t*h)
  assert abs(amp-exact)<m.mpf('1e-55')
  for p in [2,3]:
   hp=eps*m.log1p(m.mpf(1)/(p*n))/4
   overlap=m.quad(lambda y:1/m.sqrt(4*h*hp),[-hp,hp])
   assert abs(overlap-m.sqrt(hp/h))<m.mpf('1e-55')
   rows.append(dict(eps=str(eps),n=n,p=p,translation_defect=str(m.sqrt(2-2*overlap)),generator_leakage=str(h/m.sqrt(3))))
# Generic translated Gaussian orbit: bounded synthesis fails by clustering.
# Normalized Gaussian autocorrelation = exp(-a^2/4).
clusters=[]
for N in [100,10000]:
 for count in [4,16,64]:
  norm2=sum(m.exp(-m.log(m.mpf(i)/j)**2/4) for i in range(N,N+count) for j in range(N,N+count))/count
  clusters.append(dict(N=N,count=count,image_norm=str(m.sqrt(norm2)),input_norm=1))
print(json.dumps(dict(scope=__doc__,packets=rows,gaussian_clusters=clusters),indent=2))
