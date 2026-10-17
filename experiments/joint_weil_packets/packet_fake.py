"""Shifted xi-product arithmetic control, 0<theta<1/2. Non-certified."""
import mpmath as mp,json
from packet_cross import arch_cross,cross
from modulated_interval import events
mp.mp.dps=65

def fake(L,t,s,ev,theta):
    d=t-s
    if not d: prime=mp.fsum(2*mp.cosh(theta*u)*w*mp.cos(t*u) for u,w in ev)
    else: prime=mp.fsum(2*mp.cosh(theta*u)*w/(1-u/(2*L))*mp.cos((t+s)*u/2)*mp.sin(d*(L-u/2))/(L*d) for u,w in ev)
    val=-prime
    for a in [mp.mpf(1)/2+theta,mp.mpf(1)/2-theta]:
        ell=lambda v:2*mp.sinh((a+1j*v)*L)/(a+1j*v)/mp.sqrt(2*L)
        val+=arch_cross(L,t,s,a)+2*mp.re(ell(t)*ell(s))
    return val

rows=[]
for Ls in ['0.5','1.4','2']:
    L=mp.mpf(Ls); ev=events(L); nodes=[j*mp.pi/L for j in range(-8,9)]
    for t,s in [(nodes[0],nodes[0]),(nodes[1],nodes[3])]:
        assert abs(fake(L,t,s,ev,0)-2*cross(L,t,s,ev))<mp.mpf('1e-55')
    M=mp.matrix([[fake(L,t,s,ev,mp.mpf(1)/4) for s in nodes] for t in nodes])
    eigen=mp.eigsy(M,eigvals_only=True)
    rows.append(dict(L=Ls,dimension=17,theta='1/4',min_eigenvalue=mp.nstr(eigen[0],22)))
print(json.dumps(dict(scope='Shifted-product negative control, numerical only; theta=0 identity checked.',rows=rows),indent=2))
