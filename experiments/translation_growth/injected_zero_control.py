"""Positive even seed with deliberately inserted off-line transform zeros.

Uses the theta evaluation from Codex's 2026-09-30 Poisson-orbit experiment.
This is a control on a restricted translation span, not a modified zeta function.
No tabulated zeros enter the fitting. Results are floating point, not certificates.
"""
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
from scipy.linalg import lstsq, solve
from poisson_orbit_base import theta_seed

DELTA, GAMMA, TAU = .2, 20., .01
SHIFT = math.pi/GAMMA
COEFF = 1/(2*math.cosh(DELTA*SHIFT))


def seed(x, injected, dtype=np.float64):
    out = theta_seed(x, dtype)
    if injected:
        out = out + dtype(COEFF)*(theta_seed(x-dtype(SHIFT), dtype)+
                                  theta_seed(x+dtype(SHIFT), dtype))
    return out


def evaluate(centers, coeff, t, injected, dx, dtype):
    lo, hi = centers[0]-3-SHIFT, centers[-1]+3+SHIFT
    x = np.linspace(lo, hi, math.ceil((hi-lo)/dx)+1, dtype=dtype)
    c, a = np.asarray(centers,dtype=dtype), np.asarray(coeff,dtype=dtype)
    target = np.exp(-(x+dtype(t))**2/(4*dtype(TAU)))/np.sqrt(dtype(4*math.pi*TAU))
    for j in range(0,len(x),128):
        target[j:j+128] -= seed(x[j:j+128,None]-c[None,:],injected,dtype)@a
    return float(np.sqrt(simpson((target*2*np.cosh(x/2))**2,x=x)))


def lower_bound(t):
    z=np.array([complex(d,g) for d in [-DELTA,DELTA] for g in [-GAMMA,GAMMA]])
    w=z[:,None]+z.conj()[None,:]
    gram=np.ones_like(w)
    mask=abs(w)>1e-12
    gram[mask]=np.pi*w[mask]/np.sin(np.pi*w[mask])
    b=np.exp(TAU*z*z-z*t)
    return float(np.sqrt(np.vdot(b,solve(gram,b,assume_a='pos')).real))


def fit(t, injected, step):
    centers=np.linspace(-t-6,6,round((t+12)/step)+1)
    x=np.arange(centers[0]-3-SHIFT,centers[-1]+3+SHIFT,.0125)
    weight=2*np.cosh(x/2)
    target=np.exp(-(x+t)**2/(4*TAU))/math.sqrt(4*math.pi*TAU)
    mat=seed(x[:,None]-centers[None,:],injected)*weight[:,None]
    scales=np.linalg.norm(mat,axis=0)
    coeff,_,rank,singular=lstsq(mat/scales,target*weight,cond=1e-13,lapack_driver='gelsd')
    coeff/=scales
    coarse=evaluate(centers,coeff,t,injected,.0125,np.float64)
    fine=evaluate(centers,coeff,t,injected,.00625,np.longdouble)
    bound=lower_bound(t) if injected else 0.
    assert fine+1e-7>=bound,(t,fine,bound)
    assert abs(coarse-fine)<1e-5,(t,coarse,fine)
    return dict(t=t,injected=injected,step=step,residual=fine,
                grid_discrepancy=abs(coarse-fine),artificial_zero_lower_bound=bound,
                condition=float(singular[0]/singular[-1]),rank=int(rank),
                centers=centers.tolist(),coefficients=coeff.tolist())


if __name__=='__main__':
    rows=[]
    for t in [0,8,16,20,24]:
        for injected in [False,True]:
            row=fit(t,injected,1/12)
            rows.append(row)
            print(json.dumps({k:v for k,v in row.items() if k not in ['centers','coefficients']}),flush=True)
    # Refine the basis at the endpoint, preserving independently evaluated trials.
    for injected in [False,True]:
        row=fit(24,injected,1/16)
        rows.append(row)
        print(json.dumps({k:v for k,v in row.items() if k not in ['centers','coefficients']}),flush=True)
    out=dict(status='numerical control only; no RH claim',delta=DELTA,gamma=GAMMA,
             tau=TAU,shift=SHIFT,coefficient=COEFF,records=rows)
    Path(__file__).with_name('injected_zero_results.json').write_text(json.dumps(out,indent=2)+'\n')
