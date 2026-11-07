#!/usr/bin/env python3
"""Exact rational robust-positivity gate for Turbo's six-by-six Gram data.

Usage:
  python phat6_robust_interval_gate.py /path/to/phat6_data_for_gppverify.json
  python phat6_robust_interval_gate.py /path/to/data.json --margin 1/5000000 \
      --entry-error 1/100000000

The input is private GPP-bridge's MuseResearch/phat6_data_for_gppverify.json;
this public script intentionally does not copy the matrix or the research notes.

Everything asserted about the given rational Phat6 is checked over Fraction.
The conclusion concerning an UNKNOWN exact arithmetic matrix A is CONDITIONAL:
one must independently prove A is symmetric and |A_ij - Phat6_ij| <= delta.
No approximate matrix, whether rounded to 12 digits or computed from a
quadrature, comes with such an enclosure automatically.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path


def shifted_ldl_pivots(matrix: list[list[Fraction]], margin: Fraction) -> list[Fraction]:
    """Exact diagonal pivots of LDL^T of matrix - margin*I."""
    n = len(matrix)
    assert all(len(row) == n for row in matrix)
    assert all(matrix[i][j] == matrix[j][i] for i in range(n) for j in range(n))
    a = [[matrix[i][j] - (margin if i == j else 0) for j in range(n)]
         for i in range(n)]
    pivots = []
    for k in range(n):
        pivot = a[k][k]
        pivots.append(pivot)
        if pivot <= 0:
            return pivots
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] -= a[i][k] * a[k][j] / pivot
    return pivots


def check_rectangular_factor(data: dict, matrix: list[list[Fraction]]) -> None:
    if "C" not in data:
        print("Rectangular Gram factor absent from input; skipping that check")
        return
    c = [[Fraction(value) for value in row] for row in data["C"]]
    n = len(matrix)
    assert all(len(row) == n for row in c)
    for i in range(n):
        for j in range(n):
            assert sum(row[i] * row[j] for row in c) == matrix[i][j], (i, j)
    print(f"Rectangular Gram equality: {n*n}/{n*n} exact rational entries PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, help="Private bridge JSON data copied locally")
    parser.add_argument("--margin", type=Fraction, default=Fraction(1, 5_000_000))
    parser.add_argument("--entry-error", type=Fraction, default=Fraction(1, 100_000_000))
    parser.add_argument("--show-exact", action="store_true")
    args = parser.parse_args()
    if args.margin <= 0 or args.entry_error < 0:
        parser.error("margin must be positive and entry-error nonnegative")

    data = json.loads(args.data.read_text(encoding="utf-8"))
    matrix = [[Fraction(value) for value in row] for row in data["Phat6"]]
    n = len(matrix)
    assert n == 6 and all(len(row) == n for row in matrix), "expected six-by-six"
    check_rectangular_factor(data, matrix)
    pivots = shifted_ldl_pivots(matrix, args.margin)
    for k, pivot in enumerate(pivots):
        print(f"LDL pivot {k}: sign={'+' if pivot > 0 else 'NONPOSITIVE'}; {float(pivot):.15e}")
        if args.show_exact:
            print(f"  exact: {pivot}")
    if len(pivots) != n or any(p <= 0 for p in pivots):
        raise SystemExit("FAIL: chosen lower spectral margin not certified")

    # For symmetric A and E=A-Phat6, each |E_ij|<=delta implies
    # |x^T E x|<=delta*(sum_i |x_i|)^2<=n*delta*||x||^2.
    # Thus A >= (margin - n*delta)*I; this uses an INDEPENDENT error enclosure.
    lower = args.margin - n * args.entry_error
    print(f"CERTIFIED for the rational input: Phat6 >= {args.margin} * I")
    print(f"CONDITIONAL bound if arithmetic A is symmetric and |A_ij-Phat6_ij| <= "
          f"{args.entry_error}: A >= {lower} * I")
    print("ERROR ENCLOSURE IS NOT SUPPLIED OR VERIFIED BY THIS SCRIPT")
    if lower <= 0:
        raise SystemExit("Chosen entrywise error does not give a positive lower bound")


if __name__ == "__main__":
    main()
