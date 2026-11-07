"""Universal Mellin kernel for an L-fold Stieltjes dispersion representation.

If
  F(S_1,...,S_L)=int_{R_+^L} rho(x) prod_j (x_j+S_j)^(-1) d^Lx,
then, whenever Tonelli/Fubini and Mellin convergence hold and 0<Re sigma_j<1,

  M[F](sigma) = pi^L prod_j csc(pi sigma_j) M[rho](sigma).

This is the iterated identity int_0^inf S^(sigma-1)/(x+S)dS
= pi x^(sigma-1)/sin(pi sigma). It is a mathematical dispersion theorem;
it does NOT assert that an arbitrary multiloop amplitude has such an unsubtracted
multi-Stieltjes representation.

The numerical check below uses L=2 and rho(x,y)=exp(-x-y).  The one-variable
Stieltjes transform is e^S E1(S).  Its Mellin integral is evaluated directly on
[0,A], with the large-S asymptotic expansion integrated term-by-term for the tail;
this avoids overflow in e^S near the infinite endpoint.
"""
import mpmath as mp
mp.mp.dps=50


def stieltjes1(S):
    return mp.exp(S)*mp.e1(S)


def mellin_stieltjes(sig, A=mp.mpf('100'), N=45):
    head = mp.quad(lambda S: S**(sig-1)*stieltjes1(S), [0,1,10,A])
    # e^S E1(S) ~ sum_{k>=0} (-1)^k k! / S^(k+1)
    tail = mp.mpf('0')
    fact = mp.mpf(1)
    for k in range(N):
        if k:
            fact *= k
        tail += (-1)**k * fact * A**(sig-k-1)/(k+1-sig)
    return head + tail


def lhs(sig1,sig2):
    return mellin_stieltjes(sig1)*mellin_stieltjes(sig2)


def rhs(sig1,sig2):
    # Mellin rho = Gamma(sig1) Gamma(sig2).
    return ((mp.pi/mp.sin(mp.pi*sig1))*mp.gamma(sig1))* \
           ((mp.pi/mp.sin(mp.pi*sig2))*mp.gamma(sig2))


if __name__=='__main__':
    for a,b in [(mp.mpf('.3'),mp.mpf('.7')),(mp.mpf('.45'),mp.mpf('.62'))]:
        L=lhs(a,b); R=rhs(a,b)
        err=abs(L-R)/max(1,abs(R))
        assert err < mp.mpf('1e-12')
        print(a,b,'relative error',mp.nstr(err,8))
    print('PASS: L=2 Mellin transform multiplies each Stieltjes variable by pi*csc(pi*sigma)')
