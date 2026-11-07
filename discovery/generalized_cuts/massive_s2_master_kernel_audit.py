#!/usr/bin/env python3
"""Numerical audit of the exact S^2 two-affine-denominator master kernel."""

import math
import random

import mpmath as mp

mp.mp.dps = 50


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def norm2(a):
    return dot(a, a)


def qcoeff(A, B, C, D):
    q0 = C * C - norm2(D)
    dA = A - C
    dB = tuple(b - d for b, d in zip(B, D))
    q1 = 2 * (C * dA - dot(D, dB))
    q2 = dA * dA - norm2(dB)
    return mp.mpf(q0), mp.mpf(q1), mp.mpf(q2)


def reduced_parameter_integral(A, B, C, D):
    q0, q1, q2 = qcoeff(A, B, C, D)
    return 4 * mp.pi * mp.quad(lambda x: 1 / (q2 * x * x + q1 * x + q0), [0, 1])


def closed_form(A, B, C, D):
    q0, q1, q2 = qcoeff(A, B, C, D)
    disc = q1 * q1 - 4 * q2 * q0
    if abs(q2) < mp.mpf("1e-40"):
        if abs(q1) < mp.mpf("1e-40"):
            return 4 * mp.pi / q0
        return 4 * mp.pi * mp.log((q1 + q0) / q0) / q1
    if abs(disc) < mp.mpf("1e-40"):
        # Direct quadrature is the stable repeated-root limit.
        return reduced_parameter_integral(A, B, C, D)
    s = mp.sqrt(disc)
    ratio = ((2 * q2 + q1 - s) * (q1 + s)) / ((2 * q2 + q1 + s) * (q1 - s))
    return 4 * mp.pi * mp.log(ratio) / s


def direct_sphere_integral(A, B, C, D):
    def phi_integral(theta):
        st = mp.sin(theta)
        ct = mp.cos(theta)

        def f(phi):
            n = (st * mp.cos(phi), st * mp.sin(phi), ct)
            x = A + dot(B, n)
            y = C + dot(D, n)
            return st / (x * y)

        return mp.quad(f, [0, 2 * mp.pi])

    return mp.quad(phi_integral, [0, mp.pi])


def safe_case(seed):
    rng = random.Random(seed)
    B = tuple(mp.mpf(rng.uniform(-0.3, 0.3)) for _ in range(3))
    D = tuple(mp.mpf(rng.uniform(-0.3, 0.3)) for _ in range(3))
    # Choose positive offsets comfortably larger than vector norms so there are no sphere poles.
    A = mp.mpf(1.7 + rng.random())
    C = mp.mpf(1.8 + rng.random())
    return A, B, C, D


def main():
    worst_param = mp.mpf("0")
    worst_closed = mp.mpf("0")
    for seed in range(5):
        A, B, C, D = safe_case(seed)
        direct = direct_sphere_integral(A, B, C, D)
        param = reduced_parameter_integral(A, B, C, D)
        closed = closed_form(A, B, C, D)
        e1 = abs(direct - param)
        e2 = abs(param - closed)
        worst_param = max(worst_param, e1)
        worst_closed = max(worst_closed, e2)
        print(seed, mp.nstr(direct, 30), mp.nstr(e1, 5), mp.nstr(e2, 5))

    print("worst |sphere-parameter| =", mp.nstr(worst_param, 8))
    print("worst |parameter-closed| =", mp.nstr(worst_closed, 8))
    assert worst_param < mp.mpf("1e-35")
    assert worst_closed < mp.mpf("1e-35")


if __name__ == "__main__":
    main()
