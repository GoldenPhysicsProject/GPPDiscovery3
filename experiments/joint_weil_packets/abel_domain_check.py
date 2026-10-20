"""Exact rational check of the continued Abel kernel, not a claim about zeta zeros."""
import json
from fractions import Fraction as F
a=F(1,4); eps=F(1,4)
d=eps**2/(eps**2-(2*a)**2)
assert d==F(-1,3)
assert 2*(d-1)==F(-8,3)
print(json.dumps(dict(character_real_parts=['1/4','-1/4'],epsilon=str(eps),continued_diagonal=str(d),continued_off_diagonal='1',eigenvalues=[str(d+1),str(d-1)],quadratic_at_1_minus1=str(2*(d-1)),integral_convergence_requires='epsilon > 1/2',scope='Exact rational counterexample to positivity of continued Abel kernel; no zeta zero asserted'),indent=2))
