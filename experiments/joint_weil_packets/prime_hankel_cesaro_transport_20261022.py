"""Exact prime–continuum Hankel transfer from Astra's pole-balanced J_epsilon.

Mathematical identity, for eps>=0 and f compact C^2 in (0,infinity):

 H_eps f(u)=sum_(n>=2) Lambda(n)n^(-1/2-eps)f(log n-u)
              -exp((1/2-eps)u)int exp((1/2-eps)v)f(v)dv
           =int J_eps(u+v)[f''(v)-f'(v)/2]dv.

J_eps has the signed delta_at_1 boundary, so on the causal whole line:
 [d^2/dt^2 +(1/2)d/dt] J_eps(t)1_(t>=0)
 =sum_n Lambda(n)n^(-1/2-eps)delta_(log n)
  -1_(t>=0) exp((1/2-eps)t)dt -delta_0.
This is distributional; no zeta-zero / RH assumptions.
Numerical evaluation integrates at exact prime-power breakpoints.
"""
import math
from scipy.integrate import quad

def prime_powers(X):
    sieve=bytearray(X+1)
    out=[]
    for p in range(2,X+1):
        if not sieve[p]:
            if p*p<=X:
                for q in range(p*p,X+1,p): sieve[q]=1
            n=p
            while n<=X:
                out.append((math.log(n),math.log(p)))
                n*=p
    return sorted(out)

A,B=.25,.75
POWERS=prime_powers(350)

def f(v):
    if not A<v<B:return 0.
    x,y=v-A,B-v
    return x**3*y**3

def f1(v):
    if not A<v<B:return 0.
    x,y=v-A,B-v
    return 3*x*x*y*y*(y-x)

def f2(v):
    if not A<v<B:return 0.
    x,y=v-A,B-v
    return 6*x*y**3-18*x*x*y*y+6*x**3*y

def J_eps(t,eps):
    P=2*sum(lp*math.exp(-eps*ln)*(math.exp(-ln/2)-math.exp(-t/2))
            for ln,lp in POWERS if ln<=t)
    w=math.exp(-t/2)
    if abs(eps-.5)<1e-12:
        C=2*(t-1+w)
    elif abs(eps-1)<1e-12:
        C=2*(2*(1-w)-w*t+1-w)
    else:
        C=2*((math.exp((.5-eps)*t)-1)/(.5-eps)
             -w*(math.exp((1-eps)*t)-1)/(1-eps)+1-w)
    return P-C

def check(eps,u):
    direct=sum(lp*math.exp(-(.5+eps)*ln)*f(ln-u)
               for ln,lp in POWERS)
    direct-=math.exp((.5-eps)*u)*quad(
        lambda v:math.exp((.5-eps)*v)*f(v),A,B,epsabs=2e-14)[0]
    breaks=[A]+[ln-u for ln,lp in POWERS
                if A+1e-10<ln-u<B-1e-10]+[B]
    via_j=sum(quad(lambda v:J_eps(u+v,eps)*(f2(v)-.5*f1(v)),
                 x,y,epsabs=2e-14)[0]
              for x,y in zip(breaks[:-1],breaks[1:]))
    return direct,via_j,abs(direct-via_j)

if __name__=="__main__":
    errors=[]
    for eps in [0,.1,.5,1]:
        for u in [0,.5,1,2,3,4,5]:
            left,right,error=check(eps,u)
            print(f"eps={eps:<3} u={u:<3} direct={left:+.14e} J-transfer={right:+.14e} abs_error={error:.2e}")
            errors.append(error)

    # Shifted completed-xi control: exact 2cosh(theta log n) prime
    # weights equal the sum of the eps=-theta and eps=+theta transports.
    for u in [0.,1.,3.,5.]:
        theta=.25
        direct=sum(lp*math.exp(-ln/2)*2*math.cosh(theta*ln)*f(ln-u)
                   for ln,lp in POWERS)
        direct-=quad(lambda v:2*math.cosh(theta*(u+v))*
                     math.exp((u+v)/2)*f(v),
                     A,B,epsabs=2e-14)[0]
        breaks=[A]+[ln-u for ln,lp in POWERS
                     if A+1e-10<ln-u<B-1e-10]+[B]
        paired=sum(quad(lambda v:(J_eps(u+v,-theta)+J_eps(u+v,theta))*
                        (f2(v)-.5*f1(v)),x,y,epsabs=2e-14)[0]
                   for x,y in zip(breaks[:-1],breaks[1:]))
        err=abs(direct-paired)
        print(f"theta={theta} u={u} direct={direct:+.14e} paired={paired:+.14e} error={err:.2e}")
        errors.append(err)
    assert max(errors)<1e-11
    print("PASS: 32 all-prime-power comparisons, including exact F_theta 2cosh weight; max residual",max(errors))

