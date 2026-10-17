"""Actual prime-power diagnostics; 60-digit numerics, NOT interval certificates."""
import json, math
import mpmath as mp
mp.mp.dps=60
Ls=['1.4','2','3','4','5','6']
N=int(mp.ceil(mp.exp(2*mp.mpf(Ls[-1]))))
sieve=bytearray(b'\x01')*(N+1); sieve[:2]=b'\x00\x00'
for p in range(2,math.isqrt(N)+1):
    if sieve[p]: sieve[p*p:N+1:p]=b'\x00'*((N-p*p)//p+1)
events=[]
for p in range(2,N+1):
    if sieve[p]:
        n=p
        while n<=N:
            events.append((mp.log(n),mp.log(p)/mp.sqrt(n))); n*=p
events.sort()
rows=[]
for Ls0 in Ls:
    L=mp.mpf(Ls0); t=2*L
    A=mp.fsum(w*(t-u) for u,w in events if u<t)
    cp=A/L
    V=8*(mp.cosh(L)-1)
    pole=V/L
    T=V-A
    residual=cp-pole
    assert abs(residual+T/L)<mp.mpf('1e-50')
    # Autocorrelation of the normalized interval indicator is (1-u/(2L))_+.
    arch=mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi)+mp.quad(lambda u: mp.exp(-u/2)/(-mp.expm1(-2*u))*u/L,[0,2*L])+2*mp.quad(lambda u:mp.exp(-u/2)/(-mp.expm1(-2*u)),[2*L,mp.inf])
    row=dict(L=Ls0,prime_constant_rayleigh=mp.nstr(cp,18),pole=mp.nstr(pole,18),prime_minus_pole=mp.nstr(residual,18),arch=mp.nstr(arch,18),completed_constant=mp.nstr(arch-residual,18),necessary_log10_Omega=mp.nstr((mp.log(2*mp.pi)+cp)/mp.log(10),18))
    rows.append(row)
print(json.dumps({'scope':'Non-rigorous finite diagnostics. Exact identities proved separately. No RH claim.','rows':rows},indent=2))
