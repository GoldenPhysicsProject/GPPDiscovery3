"""Actual massive Euclidean double box: independent momentum/parameter QMC.
Not celestial sewing, not an on-shell scattering amplitude or rigorous bound.
All masses squared = 1. Measure d4l d4k / pi**4.
"""
import json
import numpy as np
import sympy as sp
from scipy.special import ndtri
from scipy.stats import qmc

# Symbolic successive Gaussian integration retains the shared edge.
a,b,h,v,w=sp.symbols('a b h v w', positive=True)
assert sp.simplify(a*(b-h*h/a)-(a*b-h*h)) == 0
assert sp.simplify(v*v/a+(w+h*v/a)**2/(b-h*h/a)
                   -(b*v*v+2*h*v*w+a*w*w)/(a*b-h*h)) == 0

p1=np.array([.4,.2,0.,0.]); p2=np.array([-.1,.3,.2,0.])
p4=np.array([.2,-.2,.1,.3]); P=p1+p2
def sq(x): return np.sum(x*x,axis=-1)
def parameter(z):
    # Ordered spacings give uniform Dirichlet(1,...,1), density 6!.
    z=np.sort(z,axis=1)
    x=np.diff(np.column_stack((np.zeros(len(z)),z,np.ones(len(z)))),axis=1)
    return parameter_core(x)/360

def parameter_core(x):
    h=x[:,6]; a=x[:,:3].sum(1)+h; b=x[:,3:6].sum(1)+h
    U=a*b-h*h
    v=x[:,1,None]*p1+x[:,2,None]*P
    w=-x[:,4,None]*p4+x[:,5,None]*P
    C=1+x[:,1]*sq(p1)+x[:,2]*sq(P)+x[:,4]*sq(p4)+x[:,5]*sq(P)
    F=U*C-b*sq(v)-a*sq(w)-2*h*np.sum(v*w,axis=1)
    assert np.all(F>0)
    return U/F**3

def parameter_importance(z):
    y=ndtri(z)**2
    x=y/y.sum(1)[:,None]
    # Dirichlet(1/2) density; includes Gamma(3) from Feynman integration.
    return (16*np.pi**3/15)*np.sqrt(np.prod(x,axis=1))*parameter_core(x)

def momentum(z):
    # Density 2/pi^2*(1+|l|^2)^-3 for each four-vector.
    l=ndtri(z[:,:4])/np.sqrt(-2*np.log(z[:,8,None]))
    k=ndtri(z[:,4:8])/np.sqrt(-2*np.log(z[:,9,None]))
    den=np.ones(len(z))
    for q in (l,l+p1,l+P,k,k-p4,k+P,l-k): den*=1+sq(q)
    return (1+sq(l))**3*(1+sq(k))**3/(4*den)

out={'scope':'massive off-shell Euclidean scalar double box; randomized QMC, not certified',
     'p1':p1.tolist(),'p2':p2.tolist(),'p4':p4.tolist(),'m2':1,
     'symbolic_coupled_gaussian_identity':True,'runs':[]}
for power in (18,20):
    row={'power':power,'scrambles':8}
    for name,dim,fun in [('momentum',10,momentum),('parameter',7,parameter_importance)]:
        vals=[]
        for seed in range(8):
            z=qmc.Sobol(dim,scramble=True,seed=100+seed).random_base2(power)
            vals.append(float(fun(z).mean()))
        row[name]={'mean':float(np.mean(vals)), 'scramble_se':float(np.std(vals,ddof=1)/np.sqrt(8)), 'values':vals}
    err=np.hypot(row['momentum']['scramble_se'],row['parameter']['scramble_se'])
    row['difference_over_scramble_se']=abs(row['momentum']['mean']-row['parameter']['mean'])/err
    out['runs'].append(row)
print(json.dumps(out,indent=2))
