"""Off-diagonal completed Weil forms of interval Fourier packets.
Discovery arithmetic, not a certified positivity proof.
"""
import mpmath as mp
import json
from modulated_interval import events,packet
mp.mp.dps=65

def arch_cross(L,t,s,alpha=mp.mpf(1)/2):
    L,t,s=map(mp.mpf,(L,t,s)); d=t-s; T=2*L
    alpha=mp.mpf(alpha)
    z=lambda v:alpha/2-1j*v/2
    q=mp.exp(-2*T)
    def lerch1(v):
        total=mp.mpc(0); j=0
        while True:
            total+=q**j/(j+z(v)); j+=1
            if q**j/((j+alpha/2)*(1-q))<mp.mpf('1e-60'):return total
    if not d:
        zz=z(t); j=0; phi=mp.mpc(0)
        while True:
            phi+=q**j/(j+zz)**2; j+=1
            if q**j/((j+alpha/2)**2*(1-q))<mp.mpf('1e-60'):break
        return mp.re(mp.digamma(zz))-mp.log(mp.pi)+mp.re(mp.polygamma(1,zz)-mp.exp((-alpha+1j*t)*T)*phi)/(2*T)
    H=(mp.digamma(z(s))+mp.digamma(z(-t))-2*mp.digamma(alpha/2))/2
    H+=mp.exp(-alpha*T)/2*(mp.exp(1j*s*T)*lerch1(s)+mp.exp(-1j*t*T)*lerch1(-t))
    omega=mp.digamma(alpha/2)-mp.log(mp.pi)
    return omega*mp.sin(d*L)/(d*L)+mp.im(mp.exp(1j*d*L)*H)/(L*d)

def cross(L,t,s,ev):
    L,t,s=map(mp.mpf,(L,t,s)); d=t-s
    if not d:return packet(L,t,ev)[0]
    # events stores diagonal prime weights; recover the original w_n.
    prime=mp.fsum(w/(1-u/(2*L))*mp.cos((t+s)*u/2)*mp.sin(d*(L-u/2))/(L*d) for u,w in ev)
    ell=lambda v:2*mp.sinh((mp.mpf(1)/2+1j*v)*L)/(mp.mpf(1)/2+1j*v)/mp.sqrt(2*L)
    pole=2*mp.re(ell(t)*ell(s))
    return arch_cross(L,t,s)+pole-prime

def direct_arch_cross(L,t,s,alpha=mp.mpf(1)/2):
    T=2*L; d=t-s; norm=mp.sin(d*L)/(d*L)
    k=lambda u:mp.exp(-alpha*u)/(-mp.expm1(-2*u))
    r=lambda u:2*mp.cos((t+s)*u/2)*mp.sin(d*(L-u/2))/(L*d)
    return (mp.digamma(alpha/2)-mp.log(mp.pi))*norm+mp.quad(lambda u:k(u)*(2*norm-r(u)),[T*j/100 for j in range(101)])+2*norm*mp.quad(k,[T,mp.inf])

if __name__=='__main__':
    controls=[]
    for L,t,s in [(mp.mpf('1.4'),mp.mpf(3),mp.mpf(0)),(mp.mpf(2),mp.mpf(14),mp.mpf(7)),(mp.mpf(3),mp.mpf(12),mp.mpf(-12))]:
        err=abs(arch_cross(L,t,s)-direct_arch_cross(L,t,s))
        assert err<mp.mpf('1e-55'),err
        ev=events(L)
        assert abs(cross(L,t,s,ev)-cross(L,s,t,ev))<mp.mpf('1e-55')
        controls.append(dict(L=str(L),t=str(t),s=str(s),error=mp.nstr(err,8)))
    scans=[]
    for Ls in ['1.4','2','3','4']:
        L=mp.mpf(Ls); ev=events(L); nodes=[j*mp.pi/L for j in range(-8,9)]
        mat=mp.matrix([[cross(L,t,s,ev) for s in nodes] for t in nodes])
        vals=mp.eigsy(mat,eigvals_only=True)
        scans.append(dict(L=Ls,dimension=len(nodes),min_eigenvalue=mp.nstr(vals[0],22),max_eigenvalue=mp.nstr(vals[-1,0],22)))
    print(json.dumps(dict(scope='Finite subspaces, non-certified mpmath. Not full-window positivity or RH.',controls=controls,scans=scans),indent=2))
