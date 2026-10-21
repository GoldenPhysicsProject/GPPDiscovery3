"""Numerical check of reflected Selberg identity away from critical-line zeros."""
import mpmath as m,json
m.mp.dps=60
rows=[]
for t in [1,3,10,20]:
 s=m.mpf('.5')+1j*t
 Z=m.zeta(s);Z1=m.diff(m.zeta,s);Z2=m.diff(m.zeta,s,2)
 F=-Z1/Z;Fp=F**2-Z2/Z
 A=(m.re(m.digamma(s/2))-m.log(m.pi))/2
 rhs=2*A*A-m.re(Fp)-m.re(Z2/Z)
 error=abs(abs(F)**2-rhs)
 assert abs(m.re(F)-A)<m.mpf('1e-50') and error<m.mpf('1e-50')
 rows.append({'t':t,'squared_modulus':str(abs(F)**2),'signed_selberg_term':str(m.re(Z2/Z)),'identity_error':str(error)})
print(json.dumps({'scope':__doc__,'rows':rows,'upper_bound_obtained':False},indent=2))
