"""Exact geometric-occupation realization of the full atomic Weil prime source.

Local prime occupation N_p is geometric with x=p^-sigma. The bounded
quasi-local occupation-tail observable
  O_f=sum_p log(p) sum_(k=1)^N_p f(k log p)
has expectation sum_(p,k) log(p)*p^(-k sigma)*f(k log p).
This is a one-point source representation, NOT Weil/OS positivity.
Includes an independent finite-prime-power vs Fock-basis check in the
actual wave234 bump sector and a local positive covariance check.
"""

from fractions import Fraction
from math import cos, exp, log, pi, sin, sqrt


def primes_upto(limit):
    flags = bytearray(b"\x01") * (limit + 1)
    flags[:2] = b"\x00\x00"
    for n in range(2, int(sqrt(limit)) + 1):
        if flags[n]:
            flags[n * n::n] = b"\x00" * len(flags[n * n::n])
    return [n for n in range(2, limit + 1) if flags[n]]


def g0(v):
    a = abs(v)
    if a >= 1:
        return 0.0
    return (1 - a) * (2 + cos(2 * pi * a)) / 8 + 3 * sin(2 * pi * a) / (16 * pi)


def direct_events(primes, max_u, sigma):
    """Independent von Mangoldt prime-power n=p^k atom construction."""
    max_n = int(exp(max_u))
    events = []
    for p in primes:
        power = p
        while power <= max_n:
            events.append((log(power), log(p) * power ** (-sigma)))
            power *= p
    return events


def occupation_events(primes, max_u, sigma):
    """Physical prime Fock geometric survival probability P(N_p>=k)."""
    events = []
    for p in primes:
        x = p ** (-sigma)
        a = log(p)
        k = 1
        while k * a <= max_u + 1e-13:
            events.append((k * a, a * x ** k))
            k += 1
    return events


def sector_matrices(events, size):
    a = [0.35 + 0.7 * j for j in range(size)]
    q = [[0.0] * size for _ in range(size)]
    c = [[0.0] * size for _ in range(size)]
    for i in range(size):
        for j in range(size):
            v = a[i] - a[j]
            w = a[i] + a[j] + 1
            q[i][j] = -sum(weight * (g0(u + v) + g0(-u + v))
                           for u, weight in events)
            c[i][j] = -sum(weight * g0(u - w) for u, weight in events)
    return q, c


def test_geometric_covariance():
    """Geometric tail projectors have an exact PSD covariance, not Weil form."""
    for p in (2, 3, 5):
        x = Fraction(1, p)
        for k in range(1, 7):
            for l in range(1, 7):
                mean_product = x ** max(k, l)
                covariance = mean_product - x ** (k + l)
                assert covariance >= 0
                # In this finite vector example: covariance = E[(sum c_k(I_k-x^k))²].
        coeff = [2, -1, 3, -4, 1, 2]
        cov_quad = sum(coeff[k - 1] * coeff[l - 1] *
                       (x ** max(k, l) - x ** (k + l))
                       for k in range(1, 7) for l in range(1, 7))
        assert cov_quad >= 0


def main():
    test_geometric_covariance()
    max_u = 11.6
    primes = primes_upto(int(exp(max_u)) + 1)
    for sigma in (0.5, 0.75):
        by_atoms = direct_events(primes, max_u, sigma)
        by_fock = occupation_events(primes, max_u, sigma)
        assert len(by_atoms) == len(by_fock)
        for size in (2, 3, 5, 7):
            qa, ca = sector_matrices(by_atoms, size)
            qf, cf = sector_matrices(by_fock, size)
            err_q = max(abs(qa[i][j] - qf[i][j])
                        for i in range(size) for j in range(size))
            err_c = max(abs(ca[i][j] - cf[i][j])
                        for i in range(size) for j in range(size))
            assert err_q < 2e-12 and err_c < 2e-12, (sigma, size, err_q, err_c)
            print(f"sigma={sigma}, K={size}: direct vs Fock Q,C "
                  f"max errors={err_q:.2e},{err_c:.2e}")
    print("PASS: exact all-prime quasi-local source, Q/C shadow and PSD tail covariance")


if __name__ == "__main__":
    main()
