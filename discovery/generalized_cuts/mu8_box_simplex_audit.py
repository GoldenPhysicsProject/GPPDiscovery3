from fractions import Fraction
from math import factorial


def simplex_moment(*powers: int) -> Fraction:
    """Integral over the standard simplex sum a_i=1, a_i>=0."""
    return Fraction(
        1 if not powers else __import__('functools').reduce(lambda x, y: x * y, (factorial(p) for p in powers), 1),
        factorial(len(powers) - 1 + sum(powers)),
    )


m_s2 = simplex_moment(2, 0, 2, 0)
m_t2 = simplex_moment(0, 2, 0, 2)
m_cross = simplex_moment(1, 1, 1, 1)

assert m_s2 == Fraction(1, 1260)
assert m_t2 == Fraction(1, 1260)
assert m_cross == Fraction(1, 5040)

# Integral of (s a1 a3 + t a2 a4)^2 has coefficients
# s^2*m_s2 + t^2*m_t2 + 2*s*t*m_cross.
assert m_s2 == Fraction(2, 2520)
assert 2 * m_cross == Fraction(1, 2520)

# Gamma(4-eps)/Gamma(-eps) ~ -6 eps and Gamma(-2+eps) ~ 1/(2 eps).
dimension_shift_times_pole = Fraction(-6, 2)
assert dimension_shift_times_pole == -3

assert dimension_shift_times_pole * m_s2 == Fraction(-2, 840)
assert dimension_shift_times_pole * m_t2 == Fraction(-2, 840)
assert dimension_shift_times_pole * (2 * m_cross) == Fraction(-1, 840)

print('simplex moments:', m_s2, m_t2, m_cross)
print('dimension-shift x pole:', dimension_shift_times_pole)
print('residue polynomial: -(2 s^2 + 2 t^2 + s t)/840')
