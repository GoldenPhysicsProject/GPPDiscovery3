#!/usr/bin/env python3
"""Construct zero-independent Poisson approximants to a translated Gaussian.

Each saved coefficient vector defines an actual element of the Poisson range.
The quotient-distance upper-bound formula is exact, but reported integrals
are finite numerical estimates, not interval certificates. The optional zero
table is used only afterwards for a finite dual lower bound.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
from scipy.linalg import lstsq, solve


def theta_seed(x, dtype=np.float64):
    x=np.asarray(x,dtype=dtype)
    ax=np.abs(x)
    # At |x|=3 the omitted tail is already below exp(-1200).
    u=np.exp(np.minimum(ax,dtype(3)))
    result=np.zeros_like(x)
    pi=dtype('3.141592653589793238462643383279502884')
    for n in range(1,6):
        v=n*u
        result+=(4*pi*pi*v**4-6*pi*v*v)*np.exp(-pi*v*v)
    return np.where(ax<3,np.sqrt(u)*result,0)


def residual_integral(centers,coefficients,t,tau,spacing,dtype):
    lo,hi=float(centers[0]-3),float(centers[-1]+3)
    x=np.linspace(lo,hi,int(round((hi-lo)/spacing))+1,dtype=dtype)
    c=np.asarray(centers,dtype=dtype); a=np.asarray(coefficients,dtype=dtype)
    pi=dtype('3.141592653589793238462643383279502884')
    target=np.exp(-(x+dtype(t))**2/(4*dtype(tau)))/np.sqrt(4*pi*dtype(tau))
    # Chunk to keep memory bounded, retaining extended-precision evaluation.
    residual=np.empty_like(x)
    for j in range(0,len(x),256):
        residual[j:j+256]=target[j:j+256]-theta_seed(x[j:j+256,None]-c[None,:],dtype)@a
    integrand=(residual*2*np.cosh(x/2))**2
    return float(np.sqrt(simpson(integrand,x=x)))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--zeros',type=Path)
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'research/2026-09-30_gaussian_poisson_orbit_controls.json')
    args=parser.parse_args()
    tau=.01; step=1/12; radius=6; records=[]
    for t in [0,4,8,12,16,20]:
        centers=np.linspace(-t-radius,radius,int(round((t+2*radius)/step))+1)
        x=np.arange(centers[0]-3,centers[-1]+3,.01)
        target=np.exp(-(x+t)**2/(4*tau))/math.sqrt(4*math.pi*tau)
        wroot=2*np.cosh(x/2)
        A=theta_seed(x[:,None]-centers[None,:])*wroot[:,None]
        scales=np.linalg.norm(A,axis=0)
        coeff,_,rank,singular=lstsq(A/scales,target*wroot,cond=1e-13,lapack_driver='gelsd')
        coeff/=scales
        coarse=residual_integral(centers,coeff,t,tau,.01,np.float64)
        fine=residual_integral(centers,coeff,t,tau,.005,np.longdouble)
        discrepancy=abs(coarse-fine)
        assert discrepancy<2e-5,(t,coarse,fine)
        # Analytic ambient norm, derived by Gaussian integration against
        # 4cosh^2(x/2)=2+2cosh(x).
        ambient=math.sqrt((2+2*math.exp(tau/2)*math.cosh(t))/math.sqrt(8*math.pi*tau))
        rec={'t':t,'tau':tau,'basis_count':len(centers),'numerical_rank':rank,
             'normalized_matrix_condition':float(singular[0]/singular[-1]),
             'residual_norm_float64_grid':coarse,'residual_norm_extended_precision_grid':fine,
             'grid_and_precision_discrepancy':discrepancy,'ambient_gaussian_norm':ambient,
             'max_abs_coefficient':float(np.max(np.abs(coeff))),
             'centers':centers.tolist(),'coefficients':coeff.tolist()}
        records.append(rec)
        print(json.dumps({k:v for k,v in rec.items() if k not in ['centers','coefficients']}),flush=True)
    # All coefficients above were computed before loading any zero table.
    if args.zeros:
        positive=np.loadtxt(args.zeros,max_rows=30)[:,1]
        gamma=np.r_[-positive[::-1],positive]
        delta=gamma[:,None]-gamma[None,:]
        gram=np.ones_like(delta);mask=delta!=0
        gram[mask]=np.pi*delta[mask]/np.sinh(np.pi*delta[mask])
        for rec in records:
            b=np.exp(-tau*gamma**2-1j*gamma*rec['t'])
            lower=math.sqrt(float(np.vdot(b,solve(gram,b,assume_a='pos')).real))
            rec['finite_zero_dual_lower_norm']=lower
            assert lower<rec['residual_norm_extended_precision_grid']+2e-5
    out={'status':'explicit arithmetic trial functions and finite numerical norms; no all-time bound or RH proof',
         'arithmetic_inputs':'integer lattice, Gaussian derivative phi, translations; zero table not used to construct coefficients',
         'tau':tau,'trial_function':'f(u)=sum_c a_c exp(-c/2) phi(exp(-c)u)',
         'phi':'(4*pi^2*u^4-6*pi*u^2)*exp(-pi*u^2)',
         'quotient':'L2(4*cosh(x/2)^2 dx) / closure(E(S0))',
         'zero_lower_bound':'first 30 positive tabulated ordinates and their negatives' if args.zeros else None,
         'records':records}
    args.output.write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    main()

