"""Exact rational quartet control and convergent prime-series Pick diagnostics."""
import cmath
from fractions import Fraction as F
import importlib.util
import json
import math
from pathlib import Path

spec = importlib.util.spec_from_file_location('common', Path(__file__).parents[1]/'common_spectral_variable/audit.py')
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)


def quartet(a):
    gamma = F(14)
    ys = [F(1), F(2)]
    vals = []
    for y in ys:
        re = im = F(0)
        for u in [-gamma, gamma]:
            for v in [-a, a]:
                den = (u-gamma)**2+(v-y)**2
                re += (u-gamma)/den
                im += (y-v)/den
        vals.append((re, im))
    diag = [vals[j][1]/ys[j] for j in range(2)]
    re12 = (vals[0][1]+vals[1][1])/sum(ys)
    im12 = -(vals[0][0]-vals[1][0])/sum(ys)
    det = diag[0]*diag[1]-re12**2-im12**2
    return {'shift': str(a), 'diagonal': [str(d) for d in diag],
            'det_exact': str(det), 'det_decimal': float(det)}


def run():
    positive, negative = quartet(F(0)), quartet(F(1,5))
    assert F(positive['det_exact']) > 0
    assert F(negative['det_exact']) < 0
    nmax = 100000
    terms = []
    for p in common.primes(nmax):
        n = p
        while n <= nmax:
            terms.append((math.log(n), math.log(p)))
            n *= p
    rows = []
    def logderiv(s):
        j = sum(w*cmath.exp(-s*ln) for ln,w in terms)
        return 1/s+1/(s-1)-.5*math.log(math.pi)+.5*common.psi(s/2)-j
    for height in [0,14,21,50]:
        ss = [3-1j*height, 4-1j*height]
        vals = [logderiv(s) for s in ss]
        matrix = [[(vals[i]+vals[j].conjugate())/(ss[i]+ss[j].conjugate()-1)
                   for j in range(2)] for i in range(2)]
        det = (matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]).real
        rows.append({'height':height, 'det_diagnostic':det,
                     'diagonal':[matrix[j][j].real for j in range(2)]})
    tails = {str(sigma):nmax**(1-sigma)*(math.log(nmax)/(sigma-1)+1/(sigma-1)**2)
             for sigma in [3,4]}
    result = {'control_on_axis':positive, 'control_off_axis':negative,
              'prime_cutoff':nmax, 'prime_series_tail_bounds':tails,
              'actual_arithmetic_diagnostics':rows,
              'status':'rational controls exact; arithmetic values binary64, not full interval certificates'}
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    run()
