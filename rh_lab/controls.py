"""Exact adversarial fake-zero controls for reflection/positivity claims.

Let z=s-1/2. For rational nonzero a and b, the entire polynomial
  F(s)=[z²-(a+ib)²][z²-(a-ib)²]
has real coefficients and is invariant under s -> 1-s. It has four
off-line roots at s=1/2 +/- a +/- ib. On the critical line s=1/2+it:
  F(s)=(t²-b²+a²)²+4 a² b² > 0.
Thus critical-line positivity + reflection + conjugation alone cannot
force every zero onto the critical line.
"""

import json
from fractions import Fraction as Q


def on_line(t: Q, a: Q = Q(1, 4), b: Q = Q(2)) -> Q:
    return (t*t - b*b + a*a)**2 + 4*a*a*b*b


def zero_residual(a: Q = Q(1, 4), b: Q = Q(2)) -> tuple[Q, Q]:
    # Evaluate z^4 - 2(a²-b²)z² + (a²+b²)² at z=a+ib
    zr, zi = a*a-b*b, 2*a*b
    square_r, square_i = zr*zr-zi*zi, 2*zr*zi
    return (square_r - 2*(a*a-b*b)*zr + (a*a+b*b)**2,
            square_i - 2*(a*a-b*b)*zi)


def main() -> None:
    a, b = Q(1, 4), Q(2)
    assert zero_residual(a, b) == (0, 0)
    assert 4*a*a*b*b > 0
    for t in (Q(-5), Q(-2), Q(0), Q(2), Q(5)):
        assert on_line(t, a, b) >= 4*a*a*b*b
    print(json.dumps({"control": "off-line quartet with FE symmetry and strictly positive critical-line values",
                      "off_line_re_s": ["1/4", "3/4"], "imaginary_part": "2",
                      "minimum_on_critical_line_at_t_squared": str(b*b-a*a),
                      "strict_lower_bound": str(4*a*a*b*b),
                      "exact_root_residual": ["0", "0"], "result": "PASS: weak symmetry/line positivity is insufficient"}))


if __name__ == "__main__":
    main()
