#!/usr/bin/env python3
"""Numerically cross-check the Wave234 S1 normalization and P_K parity gates.

Requires numpy/scipy. Reproduces the Wave234 source table at K=6 and checks
same-vector zeta/theta parities at K=6 and K=10. No certified full-form signs.
"""
import json

from arch_annihilator_full_s1_probe import evaluate, sector_matrix_check


def main() -> None:
    checks=[]
    for theta in (0., .25):
        v=evaluate(6,theta,500.)
        m=sector_matrix_check(6,theta,500.)
        source=({'lambda_min_plus':2.2900e-7,'lambda_min_minus':2.2145e-7}
                if theta==0 else
                {'lambda_min_plus':-3.8091e-3,'lambda_min_minus':-1.2984e-3})
        for key,reference in source.items():
            tolerance=5e-10 if theta==0 else 3e-7
            assert abs(m[key]-reference)<tolerance,(theta,key,m[key],reference)
        for parity in ('plus','minus'):
            value=v['parity'][parity]['full_over_norm']
            assert abs(value-m[f'annihilator_rayleigh_{parity}'])<1e-9
        if theta==0:
            assert abs(v['pole']['plus']-v['continuum_prime_sewing'])<1e-12
            assert all(v['parity'][p]['full_over_norm']>0 for p in ('plus','minus'))
        else:
            assert v['parity']['plus']['full_over_norm']<-.003
            assert v['parity']['minus']['full_over_norm']>+.006
        checks.append({'K':6,'theta':theta,'matrix':m,
                       'parity':v['parity'],'prime_residual':v['arithmetic_prime_residual']})

    for theta in (0., .25):
        v=evaluate(10,theta,500.)
        if theta==0:
            assert all(v['parity'][p]['full_over_norm']>5e-5
                       for p in ('plus','minus'))
            assert abs(v['pole']['plus']-v['continuum_prime_sewing'])<1e-10
        else:
            assert v['parity']['minus']['full_over_norm']<-.017
            assert v['parity']['plus']['full_over_norm']>+.020
        checks.append({'K':10,'theta':theta,'parity':v['parity'],
                       'prime_residual':v['arithmetic_prime_residual']})
    print(json.dumps({'claim':'Wave234 source-normalized numerical parity control',
                      'proof_status':'NUMERICAL ONLY; no RH proof',
                      'checks':checks},indent=2))


if __name__=='__main__':
    main()
