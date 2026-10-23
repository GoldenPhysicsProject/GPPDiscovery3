"""Ladder transfer-operator diagnostics. Quadrature is NOT an interval certificate."""
import json
import numpy as np
from fractions import Fraction as F
from numpy.polynomial.legendre import leggauss

# Exact polynomial integral after u=tanh(x/2): integral p(x)^3 dx.
logistic_cube = F(1,32)*(2-F(4,3)+F(2,5))
assert logistic_cube == F(1,30)
assert logistic_cube/3 == F(1,90)

# Check the six ordering chambers of the three points 0,x,y.
for x,y in [(3,1),(1,3),(-3,-1),(-1,-3),(2,-1),(-2,1)]:
    z=sorted([0,x,y]); a=z[1]-z[0]; b=z[2]-z[1]
    assert sorted([abs(x),abs(y),abs(x-y)])==sorted([a,b,a+b])

def P(x):
    x=np.asarray(x); z=np.pi*x
    out=np.ones_like(z,dtype=float)
    np.divide(z,np.sinh(z),out=out,where=z!=0)
    return out
rows=[]
for n in [48,80,128]:
    t,w=leggauss(n); x=6*(t+1); w=6*w
    v=np.sqrt(w*P(x)/(2*np.pi))
    T=v[:,None]*P(x[:,None]-x[None,:])*v[None,:]
    eig=np.linalg.eigvalsh(T)
    moments=[]; z=v.copy()
    for L in range(1,9):
        val=float(v@z)
        upper=1/8 if L==1 else (1/90)*(1/np.sqrt(90))**(L-2)
        lower=(1/8)*(4/45)**(L-1)
        assert lower*(1-1e-10)<=val<=upper*(1+1e-10)
        moments.append(dict(L=L,value=val,proved_lower=lower,proved_upper=upper))
        z=T@z
    rows.append(dict(nodes=n,interval=[0,12],largest_eigenvalue=float(eig[-1]),trace=float(np.trace(T)),hilbert_schmidt_squared=float(np.sum(T*T)),moments=moments))
print(json.dumps(dict(scope=__doc__,exact_M2=str(logistic_cube/3),operator_norm_upper=1/np.sqrt(90),rows=rows),indent=2))
