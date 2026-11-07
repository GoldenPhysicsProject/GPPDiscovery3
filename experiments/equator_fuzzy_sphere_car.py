#!/usr/bin/env python3
"""Exact 1D Fourier carrier -> su2 fuzzy sphere -> CAR Fock SU2-singlet census.

Combinatorics and matrices use integers/Fraction. Entropies use float logs.
This is NOT a gravitational derivation of A/(4 lP^2) or a zeta spectral model.
"""
import argparse
import json
from math import comb, log
from fractions import Fraction as Q

def matmul(A,B):
    n=len(A)
    return [[sum((A[i][k]*B[k][j] for k in range(n)),Q(0))
             for j in range(n)] for i in range(n)]

def matsum(A,B):
    return [[x+y for x,y in zip(a,b)] for a,b in zip(A,B)]

def smul(q,A):
    return [[q*x for x in a] for a in A]

def comm(A,B):
    return matsum(matmul(A,B),smul(-1,matmul(B,A)))

def equator_sl2(N):
    """N Fourier modes; [H,E]=2E, [H,F]=-2F, [E,F]=H."""
    if N < 1:
        raise ValueError("N>=1")
    n=N-1
    E=[[Q(0) for _ in range(N)] for _ in range(N)]
    F=[[Q(0) for _ in range(N)] for _ in range(N)]
    H=[[Q(0) for _ in range(N)] for _ in range(N)]
    w=[Q(1,comb(n,k)) for k in range(N)]
    for k in range(N):
        H[k][k]=n-2*k
        if k>0:E[k-1][k]=k
        if k<n:F[k+1][k]=n-k
    I=[[Q(i==j) for j in range(N)] for i in range(N)]
    assert comm(H,E)==smul(2,E)
    assert comm(H,F)==smul(-2,F)
    assert comm(E,F)==H
    j=Q(n,2)
    cas=matsum(smul(Q(1,4),matmul(H,H)),
               smul(Q(1,2),matsum(matmul(E,F),matmul(F,E))))
    assert cas==smul(j*(j+1),I)
    for k in range(1,N):
        assert k*w[k-1]==(N-k)*w[k]  # E adjoint F for positive diagonal metric
    def ev(k):return Q(n-2*k)
    # Lagrange projectors from distinct H eigenvalues, generating all matrix units.
    for k in range(N):
        for i in range(N):
            val=Q(1)
            for a in range(N):
                if a!=k:val*=(ev(i)-ev(a))/(ev(k)-ev(a))
            assert val==Q(i==k)
    return {"N":N,"spin":str(j),"metric":[str(x) for x in w],
            "Casimir":str(j*(j+1)),"operator_dimension":N*N}

def singlet_count(N):
    """dim CAR_Fock(End(V_j))^SU2 = weight0 - weight1, exact integer DP."""
    weights=[m for ell in range(N) for m in range(-ell,ell+1)]
    assert len(weights)==N*N
    d={0:1}
    for m in weights:
        nxt=d.copy()
        for weight,count in d.items():
            nxt[weight+m]=nxt.get(weight+m,0)+count
        d=nxt
    val=d[0]-d.get(1,0)
    assert 0<val<=2**(N*N)
    return val

def report(max_N=12):
    if not 1<=max_N<=24:raise ValueError("1<=max_N<=24")
    rows=[]
    for N in range(1,max_N+1):
        rep=equator_sl2(N)
        d=singlet_count(N)
        rows.append({"N":N,"spin":rep["spin"],"harmonics":sum(2*l+1 for l in range(N)),
                     "Fock_dimension":str(2**(N*N)),"singlets":str(d),
                     "singlet_fraction":d/2**(N*N),
                     "S_singlet_per_kB":log(d),
                     "entropy_density":log(d)/(N*N),
                     "fraction_times_N6":(d/2**(N*N))*N**6})
    return {
       "status":"exact finite algebra and DP; log and ratios display only",
       "Fock_max_entropy":"N^2 ln2 (nats)",
       "original_N_mode_max_entropy":"ln N (nats), NOT N^2 ln2",
       "proved_singlet_bound":"c*2^(N^2)/N^6 <= d_N <= 2^(N^2), c=2 exp(-1/24)/(3 pi^3)",
       "proved_singlet_leading_entropy":"ln d_N = N^2 ln2 + O(log N)",
       "required_gravity_matching":"a^2 = 4 pi ln2 lP^2 for a=L/N; matching is imposed, NOT derived",
       "no_go":"equator restriction misses cos(polar angle); SU2 action and CAR are extra inputs",
       "no_claim":"no Einstein dynamics, 1/4 coefficient, bulk emergence, or RH proved",
       "rows":rows,
    }

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--max-N",type=int,default=10)
    print(json.dumps(report(ap.parse_args().max_N),indent=2))
