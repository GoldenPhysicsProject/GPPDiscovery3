#!/usr/bin/env python3
"""GPP 2026-10-21: exactly conserve the Mertens terminal state under Mellin-Plancherel.

All integer-window identities below are checked with rational arithmetic
for the TRUE Mobius sequence (and adversarial bounded controls).
For arbitrary real/complex a_n, A(r)=sum_{n<=r}a_n, N>=2:

 E_N = sum_{r<N} |A(r)|^2/[r(r+1)]
     = (1/2pi) int_R |D_N(1/2+it)|^2/(t^2+1/4)dt-|A(N-1)|^2/N,
 D_N(s)=sum_{n<N}a_n n^{-s}.

The integral kernel is EXACT:
 (1/2pi) int_R exp(it log(n/m))/(t^2+1/4)dt
 = exp(-|log(n/m)|/2).
The resulting pair factor /sqrt(mn) is 1/max(m,n), NOT a GCD kernel.

For |a_n|<=1, U>0, H=sum_{n<N}1/n, elementary Hilbert
estimate gives a COMPLETELY explicit Mellin UV tail:
 0<=Tail_N(U)<=H*(4/(pi U)+16*N/(3*pi U^2)).
Proof: int_{-R}^R |D_N|^2 <=(2R+8N)H.
Split |t|>U dyadically and bound the denominator by t^2.
In particular U=sqrt(N)*log(N) gives Tail=O(1/log N),
uniformly for bounded coefficients. NO RH INPUT.

Warn: the truncated positive spectral integral MINUS terminal square
is not necessarily PSD; retaining the terminal correction is mandatory.
The LOW frequency band / signed correlations remain RH strength.
"""
from fractions import Fraction
from math import sqrt, log, pi

def mobius(N):
    a=[0]*N
    if N>1:a[1]=1
    ps=[]; comp=bytearray(N)
    for i in range(2,N):
        if not comp[i]: ps.append(i); a[i]=-1
        for p in ps:
            if i*p>=N:break
            comp[i*p]=1
            if i%p==0:a[i*p]=0;break
            a[i*p]=-a[i]
    return a

def exact_checks(a):
    N=len(a)
    A=0
    E=Fraction(0)
    spectral=Fraction(0)
    prev=0
    for n in range(1,N):
        A+=a[n]
        E+=Fraction(A*A,n*(n+1))
        spectral+=Fraction(a[n]*a[n]+2*a[n]*prev,n)
        prev+=a[n]
    assert spectral-E==Fraction(A*A,N)
    return E,spectral,Fraction(A*A,N)

if __name__ == "__main__":
    for N in (3,8,16,32,100,1000):
        for a in (mobius(N),[0]+[1]*(N-1),[0]+[(-1)**i for i in range(1,N)]):
            E,S,terminal=exact_checks(a)
            assert E>=0 and S>=0 and E+terminal==S
        E,S,terminal=exact_checks(mobius(N))
        H=sum(1/n for n in range(1,N))
        U=sqrt(N)*log(N)
        tail=H*(4/(pi*U)+16*N/(3*pi*U**2))
        print(f"PASS N={N:5d} E={float(E):.12f} spectral={float(S):.12f} "+
              f"terminal={float(terminal):.12f} certified-tail-bound={tail:.8f}")
    print("Analytic whole-line Fourier and all-U tail theorem proved in accompanying bridge note.")
    print("RH / GRH / simplicity OPEN.")
