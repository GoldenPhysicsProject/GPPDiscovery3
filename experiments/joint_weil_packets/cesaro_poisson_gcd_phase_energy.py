#!/usr/bin/env python3
"""GPP 2026-10-20: all-prime Poisson/GCD coherent-phase Cesaro energy.

For eps>0 on K_phase=product_p U(1), use independent Poisson measures
with r_p=p^(-eps). At eps=0 this is delta at identity (perfect
coherence), at eps=infinity Haar (independent random phase).
Character covariance:
  E[chi_m conjugate(chi_n)]
    = (gcd(m,n)^2/(mn))^eps.
At eps=1/2, this is the critical GCD kernel gcd(m,n)/sqrt(mn).

Let P_chi(t) = 2 sum_{p^k<=e^t} log(p)(p^(-k/2)-e^(-t/2))
                 chi_{p^k}, C(t)=2(e^(t/2)-1).
E_eps(T)=int_0^T E_mu |P_chi(t)-C(t)|^2 dt.
The script integrates t EXACTLY over prime-power cells; arithmetic
is floating point. No zeta-zero assumptions, no RH claim.

All-infinite theorem documented privately in GPP-bridge:
  * For every fixed eps>0, E_eps(T)/e^T -> 4.
  * At eps=0, E_0(T)/e^T -> 0 unconditionally by PNT.
  * At fixed T, E_eps(T)->E_0(T) as eps->0.
Thus the infinite Cesaro limit and positive Poisson regulator removal
are nonuniform. No claimed RH positivity bound.
"""
import math


def events_up_to(X):
    sieve = bytearray(X + 1)
    events = []
    for p in range(2, X + 1):
        if not sieve[p]:
            if p * p <= X:
                for j in range(p*p, X+1, p):
                    sieve[j] = 1
            k, n = 1, p
            while n <= X:
                events.append((n, p, k))
                n *= p
                k += 1
    return sorted(events)


def cell_integrals(a, b):
    # Exact integrals of 1, q, q^-1, q^2, q^-2
    # for q(t)=exp(-t/2).
    return (
        b-a,
        2*(math.exp(-a/2)-math.exp(-b/2)),
        2*(math.exp(b/2)-math.exp(a/2)),
        math.exp(-a)-math.exp(-b),
        math.exp(b)-math.exp(a),
    )


def poisson_energy(X, eps, events=None):
    events = events if events is not None else events_up_to(X)
    # E[P_chi]=2(A-qB)
    # Var(P_chi)=4(VAA-2q VAB+q² VBB), q=exp(-t/2).
    A=B=VAA=VAB=VBB=0.0
    prior={}
    prev=0.0
    pieces=[]

    def integrate(a,b):
        length,intq,intiq,intq2,intiq2=cell_integrals(a,b)
        const=A+1
        return 4*math.fsum([
            (const*const+2*B+VAA)*length,
            -2*(const*B+VAB)*intq,
            -2*const*intiq,
            (B*B+VBB)*intq2,
            intiq2,
        ])

    for n,p,k in events:
        t=math.log(n)
        pieces.append(integrate(prev,t))
        lp=math.log(p)
        r=0.0 if math.isinf(eps) else p**(-eps)
        pk=p**(-k/2)
        rk=r**k
        A+=lp*pk*rk
        B+=lp*rk

        # Covariance inside one prime:
        # cov(chi_{p^k},chi_{p^l})=r^|k-l|-r^(k+l).
        # Different primes are independent under mu_eps.
        pairs=prior.setdefault(p,[])
        ckk=1-r**(2*k)
        dAA=pk*pk*ckk
        dAB=pk*ckk
        dBB=ckk
        for l,pl in pairs:
            covariance=r**abs(k-l)-r**(k+l)
            dAA+=2*pk*pl*covariance
            dAB+=(pk+pl)*covariance
            dBB+=2*covariance
        VAA+=lp*lp*dAA
        VAB+=lp*lp*dAB
        VBB+=lp*lp*dBB
        pairs.append((k,pk))
        prev=t
    pieces.append(integrate(prev,math.log(X)))
    return math.fsum(pieces), len(events)


if __name__ == "__main__":
    for X in (1000,10000,100000):
        ev=events_up_to(X)
        print("X=",X,"prime powers=",len(ev))
        for eps in (0.0,0.001,0.01,0.1,0.5,float('inf')):
            E,n=poisson_energy(X,eps,ev)
            print(f"  eps={eps!s:<5} E={E:.12f} E/X={E/X:.9f}")
        if X==100000:
            E0,_=poisson_energy(X,0.0,ev)
            assert abs(E0-95.076289065077192)<5e-8
            Eh,_=poisson_energy(X,float('inf'),ev)
            assert abs(Eh-395451.96045338945)<5e-6
            print("  PASS: Astra exact signed energy and independent Haar cross-check")
