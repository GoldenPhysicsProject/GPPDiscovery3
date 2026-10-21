#!/usr/bin/env python3
"""GPP Oct 21: all-infinite Astra↔Turbo Möbius future/past Riesz bridge.

UNCONDITIONAL THEOREM (classical PNT at s=1):
 sum(mu(n)/n, n>=1)=0 and sum(mu(n)log(n)/n, n>=1)=-1
as ORDINARY convergent series. Consequently, for all X>=1, any real c,
 F_c(X):=sqrt(X)*sum_(n>X) mu(n)(c+log(n/X))/n
  =sqrt(X)[A(X)-c B(X)-1],
 A(X)=sum_(n<=X)mu(n)log(X/n)/n, B(X)=sum_(n<=X)mu(n)/n.
This matches the future Riesz FORM used by Turbo wave199, for the same
literal c_0; its numerical value/normalization MUST be checked by Turbo.

Plugging Astra's exact finite J representation:
 J0(log X)=2 F_c(X)+2(c-3)sqrt(X) B(X)+C(X)+2,
 C(X)=sum_(n<=X)mu(n)/sqrt(n) D(X/n).
The identity uses no RH. The future series is conditionally convergent,
not termwise absolutely! Do not truncate it without a tail certificate.

Further, D(y)=K(y)-2sqrt(y)log(y)+6sqrt(y) from Astra has the
unconditional exact limit D(infty)=-2*zeta'(1/2), and for integer N
 D(N)=-2*zeta'(1/2)-log(2*pi)/sqrt(N)
       -log(N)/(12*N**1.5)+O((1+logN)/N**3.5).
The leading constant follows by Euler-Maclaurin; direct checks below.
Neither bridge nor asymptotic proves a Möbius cancellation estimate.
"""
from math import log, sqrt
import mpmath as mp


def mobius_sieve(N):
    mu=[0]*(N+1)
    mu[1]=1
    primes=[]
    composite=bytearray(N+1)
    for i in range(2,N+1):
        if not composite[i]:
            primes.append(i);mu[i]=-1
        for p in primes:
            q=i*p
            if q>N:break
            composite[q]=1
            if i%p==0:
                mu[q]=0
                break
            mu[q]=-mu[i]
    return mu


def past_and_finite_tail(mu,X,N,c):
    S0=sum(mu[n]/n for n in range(1,N+1))
    Slog=sum(mu[n]*log(n)/n for n in range(1,N+1))
    A=sum(mu[n]*log(X/n)/n for n in range(1,X+1))
    B=sum(mu[n]/n for n in range(1,X+1))
    tail=sum(mu[n]*(c+log(n/X))/n for n in range(X+1,N+1))
    exact_finite_error=(log(X)-c)*S0-Slog-1
    assert abs((A-c*B-1)-tail-exact_finite_error)<2e-9
    return A-c*B-1,tail,exact_finite_error


def D_integer(N):
    y=mp.mpf(N)
    return (2*sum(mp.log(m)*(mp.mpf(m)**(-mp.mpf('.5'))-1/mp.sqrt(y))
                  for m in range(1,N+1))
            -2*mp.sqrt(y)*mp.log(y)+6*mp.sqrt(y))


if __name__=="__main__":
    N=1000000
    mu=mobius_sieve(N)
    for X in (10,100,1000):
        for c in (3., 4.25):
            past,tail,error=past_and_finite_tail(mu,X,N,c)
            print(f"X={X} c={c} past={past:.12g} finite-tail={tail:.12g} "
                  f"end-error={error:.9g}")
    mp.mp.dps=45
    zp=mp.diff(mp.zeta,mp.mpf('.5'))
    print("D(infinity)=-2*zeta'(1/2)=",mp.nstr(-2*zp,35))
    for n in (10,100,1000,10000):
        d=D_integer(n);q=mp.mpf(n)
        approximant=(-2*zp-mp.log(2*mp.pi)/mp.sqrt(q)
                     -mp.log(q)/(12*q**mp.mpf('1.5')))
        print(n,"D",mp.nstr(d,26),"asymptotic error",mp.nstr(d-approximant,14))
    print("PASS: exact finite past-future telescoping and D endpoint diagnostics")
