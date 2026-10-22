"""Checks of heat-trace/celestial weight dictionary. No RH or loop-amplitude claim."""
import mpmath as m,json
m.mp.dps=45
rows=[]
for v in [m.mpf('.2'),m.mpf('1'),m.mpf('4')]:
 A=m.exp(-v/2)/(-m.expm1(-2*v));B=m.exp(-3*v/2)/(-m.expm1(-2*v))
 err=abs(A+B-1/(2*m.sinh(v/2)))
 assert err<m.mpf('1e-40')
 rows.append({'v':str(v),'heat_sum_error':str(err)})
for s in [m.mpc('1.7','.4'),m.mpc('2.2','1.1')]:
 ZA=2**(-s)*m.zeta(s,m.mpf('.25'));ZB=2**(-s)*m.zeta(s,m.mpf('.75'))
 error=abs(ZA+ZB-(2**s-1)*m.zeta(s))
 assert error<m.mpf('1e-40')
 Mellin=m.quad(lambda v:v**(s-1)*(m.exp(-v/2)+m.exp(-3*v/2))/(-m.expm1(-2*v)),[0,1,m.inf])/m.gamma(s)
 assert abs(Mellin-(ZA+ZB))<m.mpf('1e-28')
 rows.append({'s':str(s),'Hurwitz_sum_error':str(error),'heat_integral_error':str(abs(Mellin-ZA-ZB))})
print(json.dumps({'scope':__doc__,'checks':rows},indent=2))
