#!/usr/bin/env python3
"""Thermal falsifier: natural SU2 angular CAR at Schwarzschild beta_H=4*pi*R.

If energy_l = l/R then beta_H energy_l=4*pi*l, and S_N=O(1), not N².
The alternate extra-soft energy_l=l/(N*R) recovers N² scaling but is
an additional dynamical assumption. No gravity or RH result derived.
"""
import argparse
import json
from math import exp, log, log1p, pi

def h_fermion(x):
    if x<0:raise ValueError("x>=0")
    y=exp(-x)
    return log1p(y)+x*y/(1+y)

def thermal_entropy(N,soft=False):
    if N<1:raise ValueError("N>=1")
    return sum((2*l+1)*h_fermion(4*pi*l/(N if soft else 1))
               for l in range(N))

def natural_uniform_bound():
    q=exp(-4*pi)
    return log(2)+q*(3-q)/(1-q)**2+4*pi*q*(3+q)/(1-q)**3

def report():
    rows=[]
    for N in (1,2,3,4,8,16,32,64,128):
        s=thermal_entropy(N)
        t=thermal_entropy(N,soft=True)
        assert s<=natural_uniform_bound()+1e-12
        rows.append({"N":N,"natural_S":s,"natural_S_over_N2":s/N**2,
                     "extra_soft_S":t,"extra_soft_S_over_N2":t/N**2})
    return {"status":"analytic uniform bound proved; numerical values for display",
       "physical_betaH":"4*pi*R (Schwarzschild, units G=ℏ=c=kB=1)",
       "natural_mode_energy":"ell/R","thermal_dimensionless_energy":"4*pi*ell",
       "analytic_uniform_ceiling":"ln 2 + q(3-q)/(1-q)^2 + 4pi q(3+q)/(1-q)^3, q=e^(-4pi)",
       "analytic_ceiling_nats":natural_uniform_bound(),
       "natural_result":"S_N bounded independent of N; no black-hole area scaling",
       "extra_soft_energy":"ell/(N R), yielding area scaling; must be derived, not assigned",
       "SU2_symmetry":"energy depends only on ell (one energy per irrep), no symmetry breaking",
       "not_claimed":"real near-horizon physics, gravity area coefficient, RH proof",
       "rows":rows}

if __name__=="__main__":
    print(json.dumps(report(),indent=2))
