"""Compare Daniel's celestial principal-series boundary proposals.

All complex arguments in this script use arithmetic s.  Celestial
Delta = 2*s.  The numerical integration tests exact Poisson identities;
it is not a proof of RH or certified interval arithmetic.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import loggamma, zeta


FIRST_ZERO_T = 14.134725141734693790457251983562


def xi(s):
    s = np.asarray(s, dtype=np.complex128)
    # xi(0)=xi(1)=1/2 by analytic continuation; evaluate these removable
    # singularities explicitly to avoid 0*infinity in the raw formula.
    singular = (s == 0) | (s == 1)
    safe_s = np.where(singular, 0.5 + 0j, s)
    value = 0.5 * safe_s * (safe_s - 1) * np.exp(
        -(safe_s / 2) * math.log(math.pi) + loggamma(safe_s / 2)
    ) * zeta(safe_s)
    return np.where(singular, 0.5 + 0j, value)


def fake_xi(s, theta=0.15):
    return xi(s + theta) * xi(s - theta)


def poisson_weights(sigma, a, b, v):
    length = b - a
    angle = math.pi * (sigma - a) / length
    numerator = math.sin(angle) / (2 * length)
    cosh = np.cosh(math.pi * v / length)
    return numerator / (cosh - math.cos(angle)), numerator / (cosh + math.cos(angle))


def simpson(y, h):
    return h / 3 * (y[0] + y[-1] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-1:2]))


def measure(func, a_delta, b_delta, t, step=0.02, cutoff=24.0):
    a = a_delta / 2
    b = b_delta / 2
    sigma = 0.5
    n = int(round(2 * cutoff / step))
    if n % 2:
        n += 1
    v = np.linspace(-cutoff, cutoff, n + 1)
    h = v[1] - v[0]
    left_kernel, right_kernel = poisson_weights(sigma, a, b, v)
    left = simpson(left_kernel * func(a + 1j * (t - v)), h)
    right = simpson(right_kernel * func(b + 1j * (t - v)), h)
    direct = complex(func(0.5 + 1j * t))
    return {
        "left": complex(left), "right": complex(right),
        "direct": direct, "residual": abs(left + right - direct),
        "conjugacy_error": abs(right - np.conj(left)),
        "mass_left": float(simpson(left_kernel, h)),
        "mass_right": float(simpson(right_kernel, h)),
    }


def main():
    for name, func in (("xi", xi), ("shifted_fake_xi", fake_xi)):
        for da, db in ((0, 2), (0.5, 1.5), (0.25, 1.5)):
            for t in (10, FIRST_ZERO_T):
                r = measure(func, da, db, t)
                print(f"{name}: Delta [{da}, {db}] t={t:.12f}")
                for k, v in r.items():
                    print(f"  {k:>16} = {v}")
    rho = 0.5 + 1j * FIRST_ZERO_T
    theta = 0.15
    print("control off-line roots from shifts:")
    for root in (rho - theta, rho + theta):
        print(f"  s={root}, |fake_xi(s)|={abs(fake_xi(root)):.3e}")


if __name__ == "__main__":
    main()
