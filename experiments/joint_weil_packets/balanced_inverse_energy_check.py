"""Finite checks of a causal energy identity; not intervals, Lean, or an RH proof."""
import json, math
from fractions import Fraction as Q
import numpy as np
from scipy.integrate import solve_ivp

# Exact rational checks of R_e(z) A_e(z)=1 and |R_e(iw)|<=1.
for e in map(Q, ['0','.01','.25','.5','1','2']):
    for z in map(Q, ['.1','1','3']):
        a=(z+e)*(z+e+Q(1,2))/(z*(z+Q(1,2)))
        r=z*(z+Q(1,2))/((z+e)*(z+e+Q(1,2)))
        assert a*r == 1
    for w in map(Q, ['0','.1','1','10']):
        assert w*w*(w*w+Q(1,4)) <= (w*w+e*e)*(w*w+(e+Q(1,2))**2)

N=1000
lam={}
for p in range(2,N+1):
    if all(p%d for d in range(2,math.isqrt(p)+1)):
        n=p
        while n<=N:
            lam[n]=math.log(p)
            n*=p

def phi(a,t):
    return t if a==0 else math.expm1(a*t)/a

rows=[]
for X in [100,1000]:
    for e in [0,.01,.25,.5,1]:
        state=np.zeros(7) # x,y,E_e,E_F,int x^2,int y^2,E_0
        W=P=W0=P0=0.
        err=0.
        points=sorted({1,X}|{n for n in lam if n<X})
        for lo,hi in zip(points,points[1:]):
            v=lam.get(lo,0.)
            W+=v*lo**(-e-.5); P+=v*lo**(-e)
            W0+=v/math.sqrt(lo); P0+=v
            def signals(t):
                decay=math.exp(-t/2)
                c=2*(phi(.5-e,t)-decay*phi(1-e,t))+2*(1-decay)
                u=2*(W-decay*P)-c
                j0=2*(W0-decay*P0)-2*math.exp(t/2)+2
                return u,j0
            def rhs(t,s):
                u,j0=signals(t)
                x,y=s[:2]
                h=u-e*x
                f=h-e*y
                return [h,h-(e+.5)*y,u*u,f*f,x*x,y*y,j0*j0]
            sol=solve_ivp(rhs,(math.log(lo),math.log(hi)),state,
                          method='DOP853',rtol=2e-12,atol=2e-13)
            assert sol.success
            state=sol.y[:,-1]
            t=math.log(hi); u,j0=signals(t)
            err=max(err,abs(u-e*state[0]-e*state[1]-math.exp(-e*t)*j0))
        x,y,Ee,EF,Ex,Ey,E0=state
        loss=e*(x*x+y*y)+e*e*Ex+e*(1+e)*Ey
        residual=abs(Ee-EF-loss)
        assert err<1e-8 and residual<1e-8
        assert E0<=math.exp(2*e*math.log(X))*Ee+1e-8
        rows.append(dict(X=X,epsilon=e,inverse_error=err,
                         dissipation_residual=residual,balanced_energy=Ee,
                         weighted_original_energy=EF,original_energy=E0,
                         nonnegative_loss=loss))
print(json.dumps(dict(scope='Exact rational multiplier samples and floating ODE checks; analytic proof is separate; no RH bound',rows=rows),indent=2))
