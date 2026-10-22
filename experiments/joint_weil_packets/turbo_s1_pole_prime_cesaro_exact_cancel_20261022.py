"""Exact Turbo wave234 S1 pole-minus-prime normalization identity.

g(x)=sin(pi*x)^2 for 0<x<1, zero else. S(r)=∫g(x)e^(r x)dx.
G0(z)=∫g(x)g(x-z)dx; h_jk(t)=G0(t-a_j-a_k-1).

For ZETA (p={1,0}) at s=a_j+a_k:
 Pole_C(s)-Prime_C(s)
   = exp(-s/2)S(-1/2)^2
       -∫ J_0(t)(h''(t)-h'(t)/2)dt.
For the shifted theta-control p={1+theta,1-theta,theta,-theta}
with prime factor 2cosh(theta log n), replace remainder by sum
r∈{-1/2+theta,-1/2-theta} exp(r*s)S(r)^2 and J_0 by J_-theta+J_theta.

This is EXACT for all bump positions s>0. No RH assertion.
Actual numerical test uses true prime powers through natural s+2 cutoff.
"""
import math
from scipy.integrate import quad

def prime_powers(X):
    sieve=bytearray(X+1);out=[]
    for p in range(2,X+1):
        if not sieve[p]:
            if p*p<=X:
                for q in range(p*p,X+1,p):sieve[q]=1
            q=p
            while q<=X:
                out.append((math.log(q),math.log(p)))
                q*=p
    return sorted(out)

def S(r):
    if abs(r)<1e-12:return .5
    return 2*math.pi**2*math.expm1(r)/(r*(r*r+4*math.pi**2))

def G0(z):
    t=abs(z)
    if t>=1:return 0.
    return (1-t)*(.25+.125*math.cos(2*math.pi*t)) +3*math.sin(2*math.pi*t)/(16*math.pi)

def G1(z):
    t=abs(z)
    if t>=1:return 0.
    v=(math.cos(2*math.pi*t)-1)/4-(math.pi/4)*(1-t)*math.sin(2*math.pi*t)
    return (-1 if z<0 else 1)*v

def G2(z):
    t=abs(z)
    if t>=1:return 0.
    return -math.pi/4*math.sin(2*math.pi*t)-math.pi**2/2*(1-t)*math.cos(2*math.pi*t)

def J(t,eps,logprimes):
    w=math.exp(-t/2)
    p=2*sum(lp*math.exp(-eps*ln)*(math.exp(-ln/2)-w)
            for ln,lp in logprimes if ln<=t)
    if abs(eps-.5)<1e-12:
        c=2*(t-1+w)
    elif abs(eps-1)<1e-12:
        c=2*(2*(1-w)-t*w+1-w)
    else:
        c=2*((math.exp((.5-eps)*t)-1)/(.5-eps)-
              w*(math.exp((1-eps)*t)-1)/(1-eps)+1-w)
    return p-c

def compare(s,theta,logprimes):
    poles=[.5,-.5] if theta==0 else [.5+theta,.5-theta,-.5+theta,-.5-theta]
    pole=sum(math.exp(r*s)*S(r)**2 for r in poles)
    weight=lambda ln: 1 if theta==0 else 2*math.cosh(theta*ln)
    prime=sum(lp*math.exp(-ln/2)*weight(ln)*G0(ln-s-1)
               for ln,lp in logprimes)
    low=[-.5] if theta==0 else [-.5+theta,-.5-theta]
    remainder=sum(math.exp(r*s)*S(r)**2 for r in low)
    edges=sorted(set([s,s+1,s+2]+[ln for ln,lp in logprimes if s<ln<s+2]))
    def filtered(t):
        return J(t,0,logprimes) if theta==0 else J(t,-theta,logprimes)+J(t,theta,logprimes)
    correction=sum(quad(lambda t:filtered(t)*(G2(t-s-1)-.5*G1(t-s-1)),
                        left,right,epsabs=1e-13)[0]
                   for left,right in zip(edges[:-1],edges[1:]))
    return pole,prime,pole-prime,remainder-correction

if __name__=="__main__":
    points=[.7,1.4,2.1,3.5]
    events=prime_powers(math.ceil(math.exp(max(points)+2))+1)
    max_err=0
    for s in points:
        for theta in [0,.25]:
            pole,prime,raw,transport=compare(s,theta,events)
            err=abs(raw-transport)
            print(f"s={s:.1f} theta={theta:.2f} pole={pole:.14e} prime={prime:.14e} raw={raw:.14e} transport={transport:.14e} error={err:.3e}")
            max_err=max(err,max_err)
    assert max_err < 3e-13
    print("PASS: direct Turbo wave234 S1 pole/prime subtraction equals causal-J transport for zeta and shifted control, max residual",max_err)
