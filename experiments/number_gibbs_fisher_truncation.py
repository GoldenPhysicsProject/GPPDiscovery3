#!/usr/bin/env python3
"""Executable truncation audit for the two-parameter number-Gibbs Fisher metric.

For
    p_n ∝ exp(-beta log n - eta (log n)^2),
with sufficient statistics X=log n and X^2, compute the normalized raw moments
m_1,...,m_4 and

    det g = (m2-m1^2)(m4-m2^2) - (m3-m1*m2)^2.

This is a numerical discovery audit only.  It deliberately distinguishes the
number-Gibbs weight from the separate von-Mangoldt Fisher weight formalized in
the arithmetic prime-gas lane.
"""

from math import exp, log


def fisher_det_truncation(N: int, beta: float, eta: float) -> float:
    raw = [0.0] * 5
    Z = 0.0
    for n in range(1, N + 1):
        x = log(n)
        w = exp(-beta * x - eta * x * x)
        Z += w
        xp = 1.0
        for r in range(5):
            raw[r] += w * xp
            xp *= x
    m = [v / Z for v in raw]
    return (m[2] - m[1] ** 2) * (m[4] - m[2] ** 2) - (m[3] - m[1] * m[2]) ** 2


def main() -> None:
    cutoffs = (10, 100, 1_000, 10_000, 100_000)
    cases = ((2.0, 0.0), (2.0, 0.15))
    for beta, eta in cases:
        print(f"beta={beta:g}, eta={eta:g}")
        previous = None
        for N in cutoffs:
            det = fisher_det_truncation(N, beta, eta)
            delta = float('nan') if previous is None else det - previous
            print(f"  N={N:>6d}  det={det:.15g}  delta={delta:.6g}")
            if not det > 0.0:
                raise AssertionError("Fisher determinant lost strict positivity")
            previous = det


if __name__ == "__main__":
    main()
