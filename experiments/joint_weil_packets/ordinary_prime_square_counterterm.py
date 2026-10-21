#!/usr/bin/env python3
"""Astra J_e = G_e + B_e + H_e: UNIFORM square-drift renormalization.

THEOREM (analytic, infinite, unconditional). Put Y=exp(t/2).
b_e(t)=t if e=0 else -expm1(-e*t)/e.
If A(Y)=sum_(p<=Y)log(p)/p obeys |A(Y)-logY|<=M,
and Chebyshev theta(Y)<=K Y, then for ALL e,t>=0,
    |B_e(t)-b_e(t)| <= 2M+2K.
Astra also proved 0<=H_e(t)<40, hence
    |J_e(t)-(G_e(t)+b_e(t))|<=2M+2K+40.
This is a uniform-in-regulator renormalization, NOT a RH proof.

Moreover, the classical second Mertens asymptotic
sum_(n<=Y)Lambda(n)/n=logY-EulerGamma+o(1)
implies for e=0 as t->infinity
 B_0(t)+H_0(t)-t ->
   Cstar= -2 EulerGamma-2 U-2+Hinf,
 U=sum_p log(p)/(p(p-1)),
 Hinf=2 sum_p sum_(k>=3) log(p) p^(-k/2).
Thus the centred oscillation can be represented by only ordinary primes
up to a known prime-independent-in-t constant plus o(1).
No claim of unconditional boundedness of G_0+t, which is RH-strength.
"""
from math import log, sqrt, expm1


def primes(N):
    s=bytearray(N+1)
    for p in range(2,int(N**0.5)+1):
        if not s[p]:
            s[p*p::p]=b'\1'*(((N-p*p)//p)+1)
    return [p for p in range(2,N+1) if not s[p]]


def b(e,t):
    return t if e==0 else -expm1(-e*t)/e


def layers(X,e,primes_list):
    t=log(X); qx=X**(-0.5); G=B=H=0.0
    for p in primes_list:
        if p>X:break
        lp=log(p); power=p; k=1
        while power<=X:
            v=2*lp*power**(-e)*(power**(-0.5)-qx)
            if k==1:G+=v
            elif k==2:B+=v
            else:H+=v
            power*=p;k+=1
    return G,B,H,b(e,t)


if __name__=="__main__":
    ps=primes(1000000)
    U=sum(log(p)/(p*(p-1)) for p in ps)
    Hinf=sum(2*log(p)*p**(-1.5)/(1-p**(-0.5)) for p in ps)
    cstar=-2*0.5772156649015328606-2*U-2+Hinf
    print("partial Cstar (NOT certified tail):",cstar)
    for X in (1000,10000,100000,1000000):
        print("X",X)
        for e in (0.,.001,.01,.1,.5):
            G,B,H,drift=layers(X,e,ps)
            print("  e",e,"B-b",round(B-drift,12),
                  "H",round(H,12),"B+H-b",round(B+H-drift,12))
        if X==100000:
            G,B,H,_=layers(X,0,ps)
            J=G+B+H-2*(sqrt(X)-1)
            assert abs(J-(-3.3544756653199615)) < 1e-8
            print("PASS: Astra J(log(1e5)) prime-only benchmark",J)
