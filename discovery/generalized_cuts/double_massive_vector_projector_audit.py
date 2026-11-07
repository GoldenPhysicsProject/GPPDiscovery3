#!/usr/bin/env python3
"""Audit the honest two-massive-vector state sum on the 5D four-gluon tree.

Codex/GPT Golden Physics discovery track.

This script sharpens the 5D Ward reconstruction result.  A raw rank-two current
J^{AB} obtained by exposing two vector legs is NOT separately Ward-transverse in one
index when the other exposed index is filled with an arbitrary basis vector.  Thus one
must not replace both massive-vector projectors by bare metrics at the unprojected
rank-two tensor level.

What is exact and safe is the physical state sum.  For the KK kinematics used here,
the three four-dimensional massive-vector polarizations can be chosen as the three
spatial unit vectors, lifted to five dimensions with epsilon^5 = 0.  We verify that

  sum_{lambda1,lambda2=1}^3 A_L(lambda1,lambda2) A_R(lambda1,lambda2)

agrees with the double four-dimensional massive-projector contraction of the exposed
rank-two tensors.  This is the correct object to use before the scalar subtraction in
C^(4) = C^(V_m) - C^(S).

Metric convention: diag(+,-,-,-,-).  All momenta are outgoing.
"""

from __future__ import annotations

import numpy as np

from massive_vector_5d_ward_reconstruction import G, amplitude, kinematics, mdot, transverse

G4 = G[:4, :4]


def exposed_two_vector_tensor(
    ks: list[np.ndarray], eps3: np.ndarray, eps4: np.ndarray
) -> np.ndarray:
    """Recover contravariant J^{AB} for exposed legs 1 and 2."""
    J = np.zeros((5, 5))
    for a in range(5):
        for b in range(5):
            e1 = np.zeros(5)
            e2 = np.zeros(5)
            e1[a] = 1.0
            e2[b] = 1.0
            # A = epsilon_{1,A} epsilon_{2,B} J^{AB}; basis lowering contributes
            # one diagonal metric sign per exposed leg.
            J[a, b] = G[a, a] * G[b, b] * amplitude(ks, [e1, e2, eps3, eps4])
    return J


def massive_projector_lower(k5: np.ndarray) -> np.ndarray:
    """Lower-index 4D massive spin-1 projector -g + p p / mu^2."""
    p = k5[:4]
    p_lower = G4 @ p
    mu2 = float(p @ G4 @ p)
    if mu2 <= 0:
        raise AssertionError("expected timelike 4D KK momentum")
    return -G4 + np.outer(p_lower, p_lower) / mu2


def rest_massive_polarizations() -> list[np.ndarray]:
    """Three physical 4D massive polarizations lifted with epsilon^5=0."""
    out: list[np.ndarray] = []
    for axis in (1, 2, 3):
        e = np.zeros(5)
        e[axis] = 1.0
        out.append(e)
    return out


def main() -> None:
    rng = np.random.default_rng(20260826)
    raw_ward_min = float("inf")
    raw_ward_max = 0.0
    projector_state_sum_max = 0.0

    for trial in range(24):
        e = 0.85 + 0.055 * trial
        theta = 0.19 + 0.037 * trial
        ks = kinematics(e, theta)
        k1, k2, k3, k4 = ks

        eps3L = transverse(k3, rng)
        eps4L = transverse(k4, rng)
        eps3R = transverse(k3, rng)
        eps4R = transverse(k4, rng)

        JL = exposed_two_vector_tensor(ks, eps3L, eps4L)
        JR = exposed_two_vector_tensor(ks, eps3R, eps4R)

        # Negative result: the raw rank-two tensor is not separately Ward transverse
        # when the other exposed leg is arbitrary/unphysical.
        ward1 = (G @ k1) @ JL
        ward2 = JL @ (G @ k2)
        raw = max(float(np.max(np.abs(ward1))), float(np.max(np.abs(ward2))))
        raw_ward_min = min(raw_ward_min, raw)
        raw_ward_max = max(raw_ward_max, raw)

        P1 = massive_projector_lower(k1)
        P2 = massive_projector_lower(k2)
        projector_sum = float(
            np.einsum("ab,ac,bd,cd->", JL[:4, :4], P1, P2, JR[:4, :4])
        )

        explicit_sum = 0.0
        for eps1 in rest_massive_polarizations():
            assert abs(mdot(k1, eps1)) < 1e-13
            for eps2 in rest_massive_polarizations():
                assert abs(mdot(k2, eps2)) < 1e-13
                AL = amplitude(ks, [eps1, eps2, eps3L, eps4L])
                AR = amplitude(ks, [eps1, eps2, eps3R, eps4R])
                explicit_sum += AL * AR

        projector_state_sum_max = max(
            projector_state_sum_max, abs(projector_sum - explicit_sum)
        )

    print(f"min raw rank-2 Ward residual: {raw_ward_min:.6e}")
    print(f"max raw rank-2 Ward residual: {raw_ward_max:.6e}")
    print(
        "max |double-projector - explicit 3x3 state sum|: "
        f"{projector_state_sum_max:.3e}"
    )

    if raw_ward_max < 1e-6:
        raise AssertionError("unexpected raw rank-two transversality")
    if projector_state_sum_max > 2e-10:
        raise AssertionError("double massive-projector state sum mismatch")

    print("PASS: raw two-vector tensor cannot be metric-contracted via naive double Ward use.")
    print("PASS: double massive projector equals the explicit nine-state physical sum.")


if __name__ == "__main__":
    main()
