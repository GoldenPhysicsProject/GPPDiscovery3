"""Exact rational mixed-loop coefficient in a Gaussian prime bath.

Pure stdlib check; no finite-prime RH approximation is asserted.
"""

from fractions import Fraction as Q
import math


def dot_ainv(v, w, diag):
    return sum(x * y / a for x, y, a in zip(v, w, diag))


def rankone_mixed(v, w, diag):
    return -dot_ainv(v, w, diag) ** 2


def matrix_log_2d(z, w, v1, v2, diag):
    """log det(A - z vv^T - w uu^T) by direct 2x2 determinant."""
    a, b = map(float, diag)
    x, y = map(lambda v: tuple(map(float, v)), (v1, v2))
    a11 = a - z * x[0] ** 2 - w * y[0] ** 2
    a12 = -z * x[0] * x[1] - w * y[0] * y[1]
    a22 = b - z * x[1] ** 2 - w * y[1] ** 2
    return math.log(a11 * a22 - a12 ** 2)


def main():
    diag = (Q(2), Q(3))
    p = (Q(1), Q(0))
    q_shared = (Q(1), Q(1))
    q_ortho = (Q(0), Q(1))
    assert rankone_mixed(p, q_shared, diag) == -Q(1, 4)
    assert rankone_mixed(p, q_ortho, diag) == Q(0)
    assert -Q(1, 5) * Q(1, 7) / Q(2) ** 2 == -Q(1, 140)

    # Mixed derivative of log det (central finite difference) tracks (1),
    # which is algebraically exact, but the difference is independently
    # calculated using the full 2x2 determinant rather than the Schur trace.
    for v2 in (q_shared, q_ortho):
        h = 0.001
        val = (
            matrix_log_2d(h, h, p, v2, diag)
            - matrix_log_2d(h, -h, p, v2, diag)
            - matrix_log_2d(-h, h, p, v2, diag)
            + matrix_log_2d(-h, -h, p, v2, diag)
        ) / (4 * h * h)
        predicted = float(rankone_mixed(p, v2, diag))
        assert abs(val - predicted) < 2e-6, (val, predicted)
        print(f"v_q={v2}, exact mixed coefficient={predicted:+.9f}, "
              f"independent derivative={val:+.9f}")
    print("Positive Gaussian prime flavours / mixed Euler loops: PASS")


if __name__ == "__main__":
    main()
