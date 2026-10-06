"""Astra: exact rational obstruction and numerical semigroup controls.

Not a Yang--Mills or RH proof. The radial profile is tested as a proposed
stationary time correlation; this does not test every celestial reconstruction.
"""
from fractions import Fraction as F
import math
import json
from pathlib import Path

t1, t2 = F(1, 4), F(1, 2)
a = 1 + (2*t1)**2
b = 1 + (t1+t2)**2
c = 1 + (2*t2)**2
excess = a*c-b*b
assert excess == F(15, 256)
# For every h>0, determinant = (a*c)^(-h) - (b*b)^(-h) < 0.
def determinant(C):
    return C(float(2*t1))*C(float(2*t2))-C(float(t1+t2))**2

results = {
    'exact_radial_product_excess': str(excess),
    'times': [str(t1), str(t2)],
    'scope': 'Scalar stationary autocorrelation / reflected Hankel test',
    'radial_hankel_determinants': {
        str(h): determinant(lambda t: (1+t*t)**(-h))
        for h in [0.1, 0.25, 0.49, 1, 2]
    },
    'positive_controls': {
        'gapless_gamma': determinant(lambda t: (1+t)**(-0.25)),
        'gapped_shifted_gamma': determinant(lambda t: math.exp(-2*t)*(1+t)**(-0.25)),
        'two_discrete_energies': determinant(lambda t: 0.3*math.exp(-t)+0.7*math.exp(-3*t)),
    },
}
assert all(v < 0 for v in results['radial_hankel_determinants'].values())
assert all(v > 0 for v in results['positive_controls'].values())
out = Path(__file__).with_name('results.json')
out.write_text(json.dumps(results, indent=2)+'\n')
print(out.read_text())
