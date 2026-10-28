"""Exact finite-interval bulk energy and celestial-cut/Carleman multiplier.

The celestial P(t)=pi*t/sinh(pi*t) is pi times the magnitude of the
off-diagonal Dirichlet-to-Neumann entry of a positive Euclidean strip
of width pi. It is also 2*D_same(t)*C_Carleman(t).
No claim is made that this alone identifies the signed Weil form.
"""

from math import cosh, pi, sinh


def strip_dtn(t):
    if t == 0:
        return (1 / pi, -1 / pi)
    x = pi * t
    return (t * cosh(x) / sinh(x), -t / sinh(x))


def celestial(t):
    return 1 if t == 0 else pi * t / sinh(pi * t)


def carleman_bulk(t):
    return pi / (2 * cosh(pi * t))


def bulk_energy(t, f0, f1, pieces=10000):
    """Midpoint integral of the exact harmonic strip interpolant."""
    dy = pi / pieces
    energy = 0.0
    for j in range(pieces):
        y = (j + 0.5) * dy
        if t == 0:
            phi = f0 * (pi - y) / pi + f1 * y / pi
            dphi = (f1 - f0) / pi
        else:
            denom = sinh(pi * t)
            phi = (f0 * sinh(t * (pi - y)) + f1 * sinh(t * y)) / denom
            dphi = t * (-f0 * cosh(t * (pi - y)) + f1 * cosh(t * y)) / denom
        energy += (dphi * dphi + t * t * phi * phi) * dy
    return energy


def main():
    for t in (0.0, 0.1, 0.4, 1.3, 3.0):
        d, off = strip_dtn(t)
        p = celestial(t)
        c = carleman_bulk(t)
        assert abs(-pi * off - p) < 1e-12
        assert abs(2 * d * c - p) < 1e-12
        assert d + off >= -1e-14 and d - off > 0
        assert abs((d * d - off * off) - (t * t)) < 1e-11
        for f0, f1 in ((1.5, -0.7), (2.0, 0.4)):
            boundary = d * (f0 * f0 + f1 * f1) + 2 * off * f0 * f1
            energy = bulk_energy(t, f0, f1)
            assert abs(energy - boundary) < 2e-6, (t, energy, boundary)
        print(f"tau={t:3.1f} strip cross={-off:.9f} "
              f"celestial P={p:.9f} Carleman bulk={c:.9f}")
    print("PASS: positive bulk strip, signed cross DtN and exact celestial P")


if __name__ == "__main__":
    main()
