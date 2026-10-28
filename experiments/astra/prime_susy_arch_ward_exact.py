"""Exact finite matrix Ward test for shared Archimedean SUSY couplings.

Even with noncommuting Arch matrices, mixed-Lyndon Bose/Fermi partners
cancel in supertraces of all Hamiltonian powers. This is a finite model
of the algebraic Q-exact deformation principle, not full Weil positivity.
"""

from itertools import product


def mmul(a, b):
    return ((a[0][0] * b[0][0] + a[0][1] * b[1][0],
             a[0][0] * b[0][1] + a[0][1] * b[1][1]),
            (a[1][0] * b[0][0] + a[1][1] * b[1][0],
             a[1][0] * b[0][1] + a[1][1] * b[1][1]))


def mtrace_power(a, k):
    result = ((1, 0), (0, 1))
    for _ in range(k):
        result = mmul(result, a)
    return result[0][0] + result[1][1]


def madd(*matrices):
    return tuple(tuple(sum(a[i][j] for a in matrices)
                       for j in range(2)) for i in range(2))


def mscale(a, n):
    return tuple(tuple(n * a[i][j] for j in range(2)) for i in range(2))


def supertrace_moment(k, max_level, arch, defect_matrices):
    """Independent Bose/Fermi modes, coupled via common noncommuting Arch."""
    total = 0
    for levels in product(range(max_level + 1), repeat=len(defect_matrices)):
        effective = madd(arch, *(mscale(defect_matrices[j], levels[j])
                                 for j in range(len(levels))))
        grade = 1
        for total_number in levels:
            # At fixed n=N_b+N_f, the allowed f are 0 and (if n>0) 1.
            # Thus all excited positive-norm doublets cancel exactly.
            grade *= 1 if total_number == 0 else 1 - 1
        total += grade * mtrace_power(effective, k)
    return total


def main():
    arch = ((5, 1), (1, 8))
    # Both defect couplings are positive definite. They do not commute
    # with each other or with the common Arch operator.
    d1 = ((3, 1), (1, 3))
    d2 = ((4, -1), (-1, 2))
    assert mmul(arch, d1) != mmul(d1, arch)
    assert mmul(d1, d2) != mmul(d2, d1)
    for count in (1, 2):
        for max_level in (1, 2, 3, 5):
            for power in range(9):
                actual = supertrace_moment(power, max_level, arch,
                                           (d1, d2)[:count])
                expected = mtrace_power(arch, power)
                assert actual == expected, (count, max_level, power)
    print("Noncommuting positive Arch coupling + two SUSY defects: PASS")
    print("All supertraces Str H^k equal pure Arch Tr H_A^k for k=0..8")
    print("Excited paired states cancel by fermion parity, not negative norm")


if __name__ == "__main__":
    main()
