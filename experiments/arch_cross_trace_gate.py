#!/usr/bin/env python3
"""Exact consequences of the six-by-six *cross* Archimedean intervals.

Given C=-sum_n L_n L_n^T, C is negative semidefinite. The JSON intervals
bound its trace and therefore norm. If an independently defined channel
B is to make B+C >= 0, necessarily tr(B)>=-tr(C). This is a conditional
compensation budget, not an assertion about any presently unknown B.
"""

import json
from fractions import Fraction
from pathlib import Path


def main() -> None:
    source = Path(__file__).with_name("arch_hankel_sine2_k6_intervals.json")
    data = json.loads(source.read_text(encoding="utf-8"))
    scale = int(data["scale"])
    mat = data["matrix"]
    n = len(mat)
    assert n == 6 and all(len(row) == n for row in mat)
    assert all(mat[i][j] == mat[j][i] for i in range(n) for j in range(n))
    lower = sum((Fraction(int(mat[i][i]["lower_scaled"]), scale)
                 for i in range(n)), Fraction(0))
    upper = sum((Fraction(int(mat[i][i]["upper_scaled"]), scale)
                 for i in range(n)), Fraction(0))
    assert lower <= upper < 0
    diagonal_compensator_lows = [Fraction(-int(mat[i][i]["upper_scaled"]), scale)
                                  for i in range(n)]
    strongest_diagonal = max(diagonal_compensator_lows)
    # Fixed rational witness close to the principal eigenvector of -C.
    # The optimization used to discover the witness is irrelevant to its
    # validity: the following inequality uses only exact integer bounds.
    witness = (1000, 692, 485, 342, 241, 170)
    witness_norm_sq = sum(x*x for x in witness)
    witness_lower = Fraction(
        sum(witness[i] * witness[j] * (-int(mat[i][j]["upper_scaled"]))
            for i in range(n) for j in range(n)),
        scale * witness_norm_sq)
    assert witness_lower > Fraction(2138, 10_000)
    # -tr(C) is necessarily >= -upper; >= each single norm contribution.
    print(json.dumps({"matrix": "Wave234 CROSS arch only", "dimension": n,
                      "trace_lower": str(lower), "trace_upper": str(upper),
                      "necessary_compensator_trace_lower": str(-upper),
                      "necessary_compensator_diagonal_lowers": [str(x) for x in diagonal_compensator_lows],
                      "rational_rayleigh_witness": list(witness),
                      "necessary_compensator_witness_rayleigh_lower": str(witness_lower),
                      "necessary_compensator_largest_eigenvalue_lower": str(max(-upper / n, strongest_diagonal, witness_lower)),
                      "conditional": "Requires B+C PSD in exactly the same normalization; no B supplied",
                      "status": "EXACT interval consequence, not full positivity or RH"}, indent=2))


if __name__ == "__main__":
    main()
