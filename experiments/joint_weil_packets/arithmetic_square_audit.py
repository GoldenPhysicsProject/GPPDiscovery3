"""Arithmetic square decomposition. Numerical diagnostics, not certified positivity."""
import contextlib, io, json
with contextlib.redirect_stdout(io.StringIO()):
    import detector_stationarity_audit as d
m=d.m
m.mp.dps=50
G=d.G
c=d.c
alpha=m.digamma(m.mpf(1)/4)-m.log(m.pi)
def audit(u,sign,R):
    u=m.mpf(u); R=m.mpf(R)
    assert R>=u+1 and int(m.exp(R))<=d.N
    def F(t):return 2*G(t)+sign*(G(t-u)+G(t+u))
    norm=F(0)
    def diff(t):return 2*(norm-F(t))
    L=u+2
    def integrand(t):
        if abs(t)<m.mpf('1e-23'):return m.mpf(0)
        return m.exp(-t/2)/(-m.expm1(-2*t))*diff(t)
    knots=sorted(set([m.mpf(0),L]+[x for x in [m.mpf(1),u-1,u,u+1] if 0<x<L]))
    r=m.exp(-L/2)
    arch=m.quad(integrand,knots)+2*norm*(m.atanh(r)+m.atan(r))
    weights=[(n,d.lam[n]/m.sqrt(n)) for n in range(2,int(m.exp(R))+1) if d.lam[n]]
    prime_squares=m.fsum(w*diff(m.log(n)) for n,w in weights)
    counterterm=2*norm*m.fsum(w for n,w in weights)
    pole=4*c*(1+sign*m.cosh(u/2))
    value=pole+alpha*norm+arch+prime_squares-counterterm
    # Independently evaluate original combined Archimedean integral.
    def original(t):
        if abs(t)<m.mpf('1e-23'): return -m.mpf('1.5')*norm
        return (2*m.exp(-2*t)*norm-2*m.exp(-t/2)*F(t))/(-m.expm1(-2*t))
    arch_original=-(m.log(m.pi)+m.euler)*norm+m.quad(original,knots)-norm*m.log1p(-m.exp(-2*L))
    direct=pole-2*m.fsum(w*F(m.log(n)) for n,w in weights)+arch_original
    assert abs(value-direct)<m.mpf('1e-38'), (value,direct)
    return {k:str(v) for k,v in dict(u=u,sign=sign,R=R,norm=norm,pole=pole,arch_squares=arch,arch_mass=alpha*norm,prime_squares=prime_squares,diagonal_subtraction=counterterm,value=value,identity_error=abs(value-direct)).items()}
rows=[audit(u,s,m.mpf(u)+r) for u in ['0.5','1','2','4'] for s in [-1,1] for r in [1,2]]
for i in range(0,len(rows),2):
    assert abs(m.mpf(rows[i]['value'])-m.mpf(rows[i+1]['value']))<m.mpf('1e-38')
print(json.dumps(dict(scope=__doc__,alpha=str(alpha),rows=rows),indent=2))
