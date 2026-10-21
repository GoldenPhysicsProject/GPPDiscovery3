"""Exact finite Mobius-energy identities and numerical scale diagnostics. Not an RH bound."""
from fractions import Fraction as F
import json, math
NMAX=1000000
mu=[1]*(NMAX+1);mu[0]=0;prime=bytearray(b'\x01')*(NMAX+1)
prime[0]=prime[1]=0
for p in range(2,NMAX+1):
    if prime[p]:
        for n in range(p,NMAX+1,p):
            mu[n]=-mu[n]
            if n>p:prime[n]=0
        for n in range(p*p,NMAX+1,p*p):mu[n]=0
checks=[]
for N in [2,3,10,37,100]:
    M=0;E=F(0)
    for n in range(1,N):
        M+=mu[n];E+=F(M*M,n*(n+1))
    D=sum((mu[n]**2*(F(1,n)-F(1,N)) for n in range(1,N)),F(0))
    C=2*sum((mu[m]*mu[n]*(F(1,n)-F(1,N)) for n in range(2,N) for m in range(1,n)),F(0))
    shift=2*sum((mu[m]*mu[m+h]*(F(1,m+h)-F(1,N)) for h in range(1,N-1) for m in range(1,N-h)),F(0))
    assert E==D+C and C==shift and C>=-D
    checks.append({'N':N,'E_exact':str(E),'identity_exact':True,'shift_identity_exact':True})
M=Q=0;E=D1=C1=A1=0.;C0=A0=0
rows=[]
for n in range(1,NMAX):
    # Cross terms pair the new coefficient with the preceding partial sum.
    C1+=2*mu[n]*M/n;C0+=2*mu[n]*M
    A1+=2*abs(mu[n])*Q/n;A0+=2*abs(mu[n])*Q
    M+=mu[n];Q+=mu[n]**2
    E+=M*M/(n*(n+1));D1+=mu[n]**2/n
    N=n+1
    if N in [100,1000,10000,100000,1000000]:
        D=D1-Q/N;C=C1-C0/N;A=A1-A0/N
        assert abs(E-D-C)<1e-9
        rows.append({'N':N,'E':E,'diagonal':D,'signed_cross':C,'absolute_pair_majorant':A,'residual':E-D-C,'E_over_logN':E/math.log(N)})
print(json.dumps({'scope':__doc__,'exact_checks':checks,'numerical_diagnostics':rows},indent=2))
