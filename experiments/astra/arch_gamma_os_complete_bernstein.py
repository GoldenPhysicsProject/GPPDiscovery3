"""Exact Archimedean gamma kinetic as a constructive OS-positive field.

For a_n = 2n + 1/2, m(t) = Re psi(1/4 + it/2) - log(pi),

  m(t) - m(0) = sum_(n>=0) 2t^2/[a_n(a_n^2+t^2)].

Finite N: A_N(q) = lam + sum_(n<N) (2/a_n) q/(q+a_n^2).
The inverse 1/A_N(q) has simple negative-real poles and POSITIVE
residues, hence an exact OS-positive Gaussian covariance on the log line.
The infinite-N form inherits OS positivity by dominated convergence.

Pure Python standard library; test command:
    python3 research/arch_gamma_os_complete_bernstein.py
"""

from __future__ import annotations

import cmath
import math


BERNOULLI = (1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66, -691 / 2730)


def digamma(z: complex, shift: int = 40) -> complex:
    """Recurrence followed by the asymptotic expansion; validation only."""
    w = z + shift
    result = cmath.log(w) - 1 / (2 * w)
    for k, b in enumerate(BERNOULLI, 1):
        result -= b / (2 * k * w ** (2 * k))
    return result - sum(1 / (z + j) for j in range(shift))


def gamma_multiplier(t: float) -> float:
    return digamma(0.25 + 0.5j * t).real - math.log(math.pi)


def gamma_tower(t: float, count: int) -> float:
    return sum(2 * t * t / (a * (a * a + t * t))
               for n in range(count) for a in (2 * n + 0.5,))


def tail_bound(t: float, count: int) -> float:
    """Rigorous integral comparison bound for the infinite remaining tower."""
    v = count + 0.25
    return t * t / 4 * (1 / v ** 3 + 1 / (2 * v * v))


def kinetic(q: float, count: int, lam: float) -> float:
    return lam + sum(2 / a * q / (q + a * a)
                     for n in range(count) for a in (2 * n + 0.5,))


def kinetic_derivative(q: float, count: int) -> float:
    return sum(2 * a / (q + a * a) ** 2
               for n in range(count) for a in (2 * n + 0.5,))


def positive_mass_spectrum(count: int, lam: float):
    """Find all poles of 1/A_N by monotonic bisection between CBF poles."""
    assert count >= 1 and lam > 0
    poles = [-((2 * n + 0.5) ** 2) for n in range(count)]
    roots = []
    for j in range(count):
        # Left side of the interval is immediately to the RIGHT of pole j.
        # Right side is immediately to the LEFT of pole j-1, or zero.
        left = poles[j] + 1e-10
        right = poles[j - 1] - 1e-10 if j else 0.0
        assert kinetic(left, count, lam) < 0
        assert kinetic(right, count, lam) > 0
        for _ in range(110):
            middle = (left + right) / 2
            if kinetic(middle, count, lam) < 0:
                left = middle
            else:
                right = middle
        root = (left + right) / 2
        mass = math.sqrt(-root)
        residue = 1 / kinetic_derivative(root, count)
        assert residue > 0
        roots.append((mass, residue))
    contact = 1 / (lam + sum(2 / (2 * n + 0.5) for n in range(count)))
    return contact, roots


def cov_positive_time(t: float, modes) -> float:
    """Non-contact part C_N(t) for t>0."""
    assert t > 0
    return sum(residue / (2 * mass) * math.exp(-mass * t)
               for mass, residue in modes)


def os_quadratic(times: list[float], values: list[complex], modes) -> float:
    assert len(times) == len(values) and min(times) > 0
    return sum(
        residue / (2 * mass)
        * abs(sum(c * math.exp(-mass * u) for u, c in zip(times, values))) ** 2
        for mass, residue in modes
    )


def main() -> None:
    m0 = gamma_multiplier(0.0)
    assert abs(m0 - (-5.372183419225665)) < 2e-13
    print(f"Arch m(0)={m0:.15f}")
    for t in (0.1, 0.5, 1.0, 5.0, 17.0):
        target = gamma_multiplier(t) - m0
        for count in (8, 64, 1024):
            approx = gamma_tower(t, count)
            err = target - approx
            bound = tail_bound(t, count)
            assert -2e-12 < err < bound + 2e-12, (t, count, err, bound)
        print(f"t={t:5}: gamma kinetic={target:.12f}, "
              f"N1024 error={target - gamma_tower(t, 1024):+.3e}")

    for count in (1, 2, 4, 8, 16):
        lam = 1.0
        contact, modes = positive_mass_spectrum(count, lam)
        for q in (0.0, 0.1, 1.0, 9.0, 123.0):
            exact = 1 / kinetic(q, count, lam)
            reconstructed = contact + sum(res / (q + mass ** 2)
                                          for mass, res in modes)
            assert abs(exact - reconstructed) < 2e-12, (count, q)
        times = [0.12, 0.6, 1.1, 2.0]
        values = [1.3 + 0.7j, -1.0 + 1.2j, 0.3 - 0.9j, -0.4 + 0.1j]
        square = os_quadratic(times, values, modes)
        direct = sum(values[j] * values[k].conjugate()
                     * cov_positive_time(times[j] + times[k], modes)
                     for j in range(len(times)) for k in range(len(times)))
        assert abs(square - direct.real) < 1e-12
        assert abs(direct.imag) < 1e-12
        assert square > 0
        print(f"N={count:2}: contact={contact:.9f}, "
              f"m_min={min(m for m, _ in modes):.8f}, "
              f"res_min={min(w for _, w in modes):.3e}, "
              f"OS square={square:.12f}")
    print("Complete-Bernstein arch tower / positive OS residues: PASS")


if __name__ == "__main__":
    main()
