#!/usr/bin/env python3
"""Numerical audit of massive-vector polarization completeness and sewing.

Mostly-minus metric g = diag(1,-1,-1,-1).  Starting from the three rest-frame
polarizations, boost them to a generic on-shell momentum p^2 = m^2 and verify

    sum_lambda eps^mu_lambda eps^nu_lambda
      = -g^{mu nu} + p^mu p^nu / m^2.

Then generate generic currents and project them transverse to p.  The full
three-polarization state sum is checked against the exact conserved-current
collapse -J_L . J_R.
"""

from __future__ import annotations

import math
import random


def mdot(a, b):
    return a[0] * b[0] - sum(a[i] * b[i] for i in range(1, 4))


def lower(a):
    return [a[0], -a[1], -a[2], -a[3]]


def boost_polarizations(p, mass):
    energy = p[0]
    v = [p[i] / energy for i in range(1, 4)]
    v2 = sum(x * x for x in v)
    gamma = energy / mass
    pols = []

    for a in range(3):
        eps = [0.0] * 4
        eps[0] = gamma * v[a]
        for i in range(3):
            correction = 0.0 if v2 == 0.0 else (gamma - 1.0) * v[i] * v[a] / v2
            eps[i + 1] = (1.0 if i == a else 0.0) + correction
        pols.append(eps)
    return pols


def completeness_tensor(pols):
    return [
        [sum(eps[mu] * eps[nu] for eps in pols) for nu in range(4)]
        for mu in range(4)
    ]


def target_projector(p, mass):
    gdiag = [1.0, -1.0, -1.0, -1.0]
    return [
        [
            (-(gdiag[mu]) if mu == nu else 0.0) + p[mu] * p[nu] / mass**2
            for nu in range(4)
        ]
        for mu in range(4)
    ]


def transverse_part(j, p, mass):
    coeff = mdot(p, j) / mass**2
    return [j[mu] - coeff * p[mu] for mu in range(4)]


def contract_cov_projector_cov(j_left, projector, j_right):
    jl = lower(j_left)
    jr = lower(j_right)
    return sum(
        jl[mu] * projector[mu][nu] * jr[nu]
        for mu in range(4)
        for nu in range(4)
    )


def run(seed=20260825, trials=100, tol=1e-11):
    rng = random.Random(seed)
    worst = {
        "onshell": 0.0,
        "transversality": 0.0,
        "orthonormality": 0.0,
        "completeness": 0.0,
        "sewing": 0.0,
    }

    for _ in range(trials):
        mass = rng.uniform(0.3, 3.0)
        spatial = [rng.uniform(-4.0, 4.0) for _ in range(3)]
        energy = math.sqrt(mass**2 + sum(x * x for x in spatial))
        p = [energy, *spatial]

        pols = boost_polarizations(p, mass)
        projector = completeness_tensor(pols)
        target = target_projector(p, mass)

        worst["onshell"] = max(worst["onshell"], abs(mdot(p, p) - mass**2))
        worst["transversality"] = max(
            worst["transversality"],
            max(abs(mdot(p, eps)) for eps in pols),
        )
        worst["orthonormality"] = max(
            worst["orthonormality"],
            max(
                abs(mdot(pols[i], pols[j]) + (1.0 if i == j else 0.0))
                for i in range(3)
                for j in range(3)
            ),
        )
        worst["completeness"] = max(
            worst["completeness"],
            max(abs(projector[mu][nu] - target[mu][nu]) for mu in range(4) for nu in range(4)),
        )

        raw_left = [rng.uniform(-2.0, 2.0) for _ in range(4)]
        raw_right = [rng.uniform(-2.0, 2.0) for _ in range(4)]
        j_left = transverse_part(raw_left, p, mass)
        j_right = transverse_part(raw_right, p, mass)

        sewn = contract_cov_projector_cov(j_left, projector, j_right)
        collapsed = -mdot(j_left, j_right)
        worst["sewing"] = max(worst["sewing"], abs(sewn - collapsed))

    print("massive-vector projector audit (mostly-minus metric)")
    for key, value in worst.items():
        print(f"  worst {key:16s}: {value:.3e}")

    if any(value > tol for value in worst.values()):
        raise SystemExit("FAIL: tolerance exceeded")
    print(f"PASS: all {trials} trials below tolerance {tol:.1e}")


if __name__ == "__main__":
    run()
