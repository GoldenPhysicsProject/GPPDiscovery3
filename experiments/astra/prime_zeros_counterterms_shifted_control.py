"""Independent RH control: Riesz-smoothed explicit formula with OFF-LINE zeros.

For G_theta(s)=zeta(s-theta)*zeta(s+theta), compute the prime-power
Riesz sum directly and compare it to pole/zero/trivial/kernel residues.
The zeros of G_theta are rho+theta and rho-theta, off the 1/2 line.
Only the first N zeta zeros are used; the Riesz kernel decays as |rho|^-3.
"""

from __future__ import annotations

import math

import mpmath as mp


def primes_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (
                (limit - p * p) // p + 1
            )
    return [i for i in range(2, limit + 1) if sieve[i]]


def prime_riesz(x: float, a: float) -> float:
    """Sum_{p^k<x} log(p) * (p^k)^a * (1-p^k/x)^2."""
    ans = 0.0
    for p in primes_to(math.floor(x)):
        n = p
        while n < x:
            ans += math.log(p) * (n**a) * ((1 - n / x) ** 2)
            n *= p
    return ans


def log_derivative(w: mp.mpc) -> mp.mpc:
    return -mp.diff(mp.zeta, w) / mp.zeta(w)


def riesz_residue(x: mp.mpf, z: mp.mpc) -> mp.mpc:
    return 2 * (x**z) / (z * (z + 1) * (z + 2))


def explicit_components(x: float, a: float, zeros: list[mp.mpc], nzeros: int):
    """Pole, zero and small terms of the exact m=2 Riesz explicit formula."""
    x = mp.mpf(str(x))
    a = mp.mpf(str(a))
    pole = riesz_residue(x, 1 + a)
    zero_sum = -2 * sum(
        (riesz_residue(x, rho + a)).real for rho in zeros[:nzeros]
    )
    trivial = -sum(riesz_residue(x, a - 2 * k) for k in range(1, 12))
    small = (
        log_derivative(-a)
        - 2 * x ** (-1) * log_derivative(-1 - a)
        + x ** (-2) * log_derivative(-2 - a)
    )
    return tuple(float(mp.re(y)) for y in (pole, zero_sum, trivial, small))


def main() -> None:
    mp.mp.dps = 35
    theta = 0.15
    zeros = [mp.zetazero(k) for k in range(1, 161)]
    for x in (1000.3, 10000.3, 100000.3):
        prime_total = prime_riesz(x, theta) + prime_riesz(x, -theta)
        print(f"\nx={x}, theta={theta}; G_theta zeros at Re(s)=0.35 and 0.65")
        print(f"literal two-species prime-power Riesz total = {prime_total:.10f}")
        for n in (20, 40, 80, 160):
            pos = explicit_components(x, theta, zeros, n)
            neg = explicit_components(x, -theta, zeros, n)
            pole, zsum, triv, small = [p + q for p, q in zip(pos, neg)]
            predicted = pole + zsum + triv + small
            print(
                f"Nzero={n:3} pole={pole:+.7f} zero={zsum:+.7f}"
                f" trivial={triv:+.7f} small={small:+.7f}"
                f" predicted={predicted:+.7f}"
                f" residual={predicted-prime_total:+.7f}"
            )


if __name__ == "__main__":
    main()
