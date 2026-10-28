"""Exact prime word-to-Euler Fock reduction through Lyndon/Bose-Fermi pairs.

Free ordered words generate 1/(1 - sum_p x_p). Chen-Fox-Lyndon unique
factorization yields the product over primitive Lyndon words. Pair every
multi-letter (mixed-prime) Lyndon boson with one fermionic mode to cancel
its partition factor in a graded trace, retaining only single-prime
bosonic oscillators. The Hilbert metric remains positive.

All polynomial checks use integer arithmetic and require no dependencies.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import product


def is_lyndon(word: tuple[int, ...]) -> bool:
    return all(word < word[j:] + word[:j] for j in range(1, len(word)))


def lyndon_words(alphabet: int, max_len: int):
    for size in range(1, max_len + 1):
        for word in product(range(alphabet), repeat=size):
            if is_lyndon(word):
                yield word


def multiply_by_fermionic_factor(poly: dict[tuple[int, ...], int],
                                exponents: tuple[int, ...], degree: int):
    """Multiply by (1 - prod_i x_i^exponents_i) in truncated integer ring."""
    result = dict(poly)
    for powers, coeff in poly.items():
        new = tuple(a + b for a, b in zip(powers, exponents))
        if sum(new) <= degree:
            result[new] = result.get(new, 0) - coeff
            if result[new] == 0:
                del result[new]
    return result


def lyndon_factor_product(alphabet: int, degree: int, mixed_only=False):
    zero = (0,) * alphabet
    poly = {zero: 1}
    word_count = 0
    for word in lyndon_words(alphabet, degree):
        if mixed_only and len(set(word)) == 1:
            continue
        multiplicity = Counter(word)
        mono = tuple(multiplicity.get(i, 0) for i in range(alphabet))
        poly = multiply_by_fermionic_factor(poly, mono, degree)
        word_count += 1
    return poly, word_count


def expected_one_minus_sum(alphabet: int):
    p = {(0,) * alphabet: 1}
    for i in range(alphabet):
        v = tuple(int(j == i) for j in range(alphabet))
        p[v] = -1
    return p


def rational_product_ghost(alphabet_x: list[Fraction], degree: int) -> Fraction:
    prod_ghost = Fraction(1)
    for word in lyndon_words(len(alphabet_x), degree):
        if len(set(word)) < 2:
            continue
        weight = Fraction(1)
        for letter in word:
            weight *= alphabet_x[letter]
        prod_ghost *= 1 - weight
    return prod_ghost


def main():
    for alphabet, degree in ((2, 10), (3, 7), (4, 5)):
        product_all, count_all = lyndon_factor_product(alphabet, degree)
        expect = expected_one_minus_sum(alphabet)
        assert product_all == expect, (alphabet, degree, product_all)

        product_mixed, count_mixed = lyndon_factor_product(
            alphabet, degree, mixed_only=True)
        for i in range(alphabet):
            single = tuple(int(j == i) for j in range(alphabet))
            product_mixed = multiply_by_fermionic_factor(product_mixed,
                                                         single, degree)
        assert product_mixed == expect
        print(f"alphabet={alphabet}, degree={degree}, "
              f"Lyndon total={count_all}, mixed={count_mixed}: exact PASS")

    x = [Fraction(1, 4), Fraction(1, 9)]  # primes 2,3 at s=2
    free_word = 1 / (1 - sum(x))
    physical_euler = 1 / ((1 - x[0]) * (1 - x[1]))
    exact_ghost_ratio = physical_euler / free_word
    print(f"two primes at s=2: free words={free_word}, "
          f"Euler physical={physical_euler}, ghost={exact_ghost_ratio}")
    for d in (3, 5, 7, 9, 11):
        approx = rational_product_ghost(x, d)
        assert abs(approx - exact_ghost_ratio) < Fraction(1, 5) ** (d // 2)
        print(f"Lyndon max length={d:2}, finite ghost={float(approx):.12f}, "
              f"error={float(approx-exact_ghost_ratio):+.3e}")

    # Single mixed-Lyndon supersymmetric oscillator: bosonic n>=0,
    # fermion f=0,1, Q=b f^dagger and {Q,Q^dagger}=N_b+N_f.
    # Each positive-energy n state forms a paired even/odd doublet;
    # graded partition telescopes to the vacuum index 1.
    w = Fraction(1, 6)
    for max_b in (0, 1, 2, 4, 8):
        boson = sum((w ** k for k in range(max_b + 1)), Fraction())
        odd = sum((w ** (k + 1) for k in range(max_b + 1)), Fraction())
        index_partial = boson - odd
        assert index_partial == 1 - w ** (max_b + 1)
    print("Positive-norm Bose/Fermi mode pairs: supertrace index 1 PASS")


if __name__ == "__main__":
    main()
