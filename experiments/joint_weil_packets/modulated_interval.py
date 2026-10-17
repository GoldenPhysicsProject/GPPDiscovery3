"""Joint prime/pole/arch test on normalized modulated intervals.
Analytic packet formula; mpmath values are NOT interval certificates or RH.
"""
import math, json
import mpmath as mp
mp.mp.dps=45

def events(L):
    N=int(mp.ceil(mp.exp(2*L))); sieve=bytearray(b'\1')*(N+1)
    sieve[:2]=b'\0\0'
    for p in range(2,math.isqrt(N)+1):
        if sieve[p]: sieve[p*p:N+1:p]=b'\0'*((N-p*p)//p+1)
    out=[]
    for p in range(2,N+1):
        if sieve[p]:
            n=p
            while n<N:
                u=mp.log(n)
                if u<2*L: out.append((u,2*mp.log(p)/mp.sqrt(n)*(1-u/(2*L))))
                n*=p
    return out

def packet(L,tau,ev):
    L=mp.mpf(L); tau=mp.mpf(tau); T=2*L
    z=mp.mpf(1)/4-1j*tau/2
    # Exponentially convergent Lerch series. Bound its omitted absolute tail.
    q=mp.exp(-2*T); J=1
    while q**J/(J+mp.mpf(1)/4)**2/(1-q)>mp.power(10,-mp.mp.dps+5): J+=1
    lerch=mp.fsum(q**j/(j+z)**2 for j in range(J))
    correction=mp.re(mp.polygamma(1,z)-mp.exp((-mp.mpf(1)/2+1j*tau)*T)*lerch)/(2*T)
    arch=mp.re(mp.digamma(z))-mp.log(mp.pi)+correction
    pole=4/L*mp.re(mp.sinh((mp.mpf(1)/2+1j*tau)*L)**2/(mp.mpf(1)/2+1j*tau)**2)
    prime=mp.fsum(w*mp.cos(tau*u) for u,w in ev)
    return arch+pole-prime, arch, pole, prime

def direct_arch(L,tau):
    T=2*L
    k=lambda u:mp.exp(-u/2)/(-mp.expm1(-2*u))
    g=lambda u:k(u)*2*(1-(1-u/T)*mp.cos(tau*u))
    pts=[T*j/100 for j in range(101)]
    return mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi)+mp.quad(g,pts)+2*mp.quad(k,[T,mp.inf])

if __name__=='__main__':
    controls=[]
    for Ls,ts in [('1.4','0'),('1.4','3'),('2','14'),('3','30')]:
        L=mp.mpf(Ls); tau=mp.mpf(ts); ev=events(L)
        vals=packet(L,tau,ev)
        error=abs(vals[1]-direct_arch(L,tau))
        assert error<mp.mpf('1e-35'), error
        # Independently integrate the polar moment; the full pole is 2 Re(ell_+^2).
        ell=mp.quad(lambda x:mp.exp((mp.mpf(1)/2+1j*tau)*x),[-L,L])/mp.sqrt(2*L)
        assert abs(vals[2]-2*mp.re(ell**2))<mp.mpf('1e-35')
        controls.append(dict(L=Ls,tau=ts,arch_error=mp.nstr(error,8)))
    scans=[]
    for Ls in ['1.4','2','3','4']:
        L=mp.mpf(Ls); ev=events(L)
        grid=[mp.mpf(j)/2 for j in range(81)]+[mp.power(10,mp.mpf(j)/8) for j in range(13,49)]
        best=min((packet(L,t,ev)[0],t) for t in grid)
        scans.append(dict(L=Ls,tests=len(grid),min_Q=mp.nstr(best[0],18),tau=mp.nstr(best[1],18)))
    print(json.dumps(dict(scope='Finite packet sampling only; no arbitrary-test positivity or RH.',controls=controls,scans=scans),indent=2))
