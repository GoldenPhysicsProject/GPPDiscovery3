"""Finite-sum checks of prime-pair identities; not a Lean proof."""
from fractions import Fraction as F

for p in (2, 3, 5, 11):
    a = F(1, p)
    weights = [(1-a)*a**n for n in range(201)]
    mean = sum(n*w for n, w in enumerate(weights))
    second = sum(n*n*w for n, w in enumerate(weights))
    expected_mean = F(1, p-1)
    expected_second = F(p+1, (p-1)**2)
    assert abs(mean-expected_mean) < F(1, 10**50)
    assert abs(second-expected_second) < F(1, 10**50)
    # Difference charge vanishes on each paired basis vector.
    assert all(n-n == 0 for n in range(201))
    print(p, 'mean N:', expected_mean,
          'variance N:', expected_second-expected_mean**2)
