#!/usr/bin/env python3
"""Audit the correct Ward structure of the 5D -> 4D massive-vector cut state sum.

Codex/GPT Golden Physics discovery track.

The earlier conditional shortcut p_4D . J_4D = 0 is generically false for a massive
vector obtained from a five-dimensional massless gauge field.  The exact identity is

    K . J_5D = 0,  K = (p, kappa)
        => p . J_4D = kappa J_5.

Consequently, with mu^2 = kappa^2,

    -J_L,4 . J_R,4 + (p.J_L,4)(p.J_R,4)/mu^2
      = -J_L,4 . J_R,4 + J_L,5 J_R,5
      = -J_L,5D . J_R,5D.

This script checks the statement on the complete color-ordered 5D four-gluon tree
built from the standard cubic and quartic Yang-Mills vertices.  It verifies full 5D
Ward identities while exhibiting nonzero 4D longitudinal contractions.

Metric convention: diag(+,-,-,-,-). All external momenta are outgoing.
"""

from __future__ import annotations

import numpy as np

G = np.diag([1.0, -1.0, -1.0, -1.0, -1.0])


def mdot(a: np.ndarray, b: np.ndarray) -> float:
    return float(a @ G @ b)


def v3(p: np.ndarray, q: np.ndarray, r: np.ndarray) -> np.ndarray:
    """Color-ordered cubic vertex V^{mu nu rho}(p,q,r), all momenta outgoing."""
    d = len(p)
    out = np.zeros((d, d, d))
    for mu in range(d):
        for nu in range(d):
            for rho in range(d):
                out[mu, nu, rho] = (
                    G[mu, nu] * (p - q)[rho]
                    + G[nu, rho] * (q - r)[mu]
                    + G[rho, mu] * (r - p)[nu]
                )
    return out


def transverse(k: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Generate a generic vector transverse to the null momentum k."""
    v = rng.normal(size=5)
    n = np.array([1.0, 0.0, 0.0, 0.0, 0.0])
    denom = mdot(k, n)
    if abs(denom) < 1e-14:
        raise RuntimeError("reference vector accidentally orthogonal to k")
    return v - n * (mdot(k, v) / denom)


def four_gluon_terms(ks: list[np.ndarray], eps: list[np.ndarray]) -> tuple[float, ...]:
    """Return the two exchange terms and three quartic metric structures.

    External polarizations are contravariant.  Since V3 is written with upper
    Lorentz indices, external polarizations are lowered before contraction.
    """
    k1, k2, k3, k4 = ks
    e1, e2, e3, e4 = [G @ e for e in eps]

    p12 = k1 + k2
    a12 = np.einsum(
        "m,n,mna,ab,brs,r,s->",
        e1,
        e2,
        v3(k1, k2, -p12),
        G,
        v3(p12, k3, k4),
        e3,
        e4,
    ) / mdot(p12, p12)

    p23 = k2 + k3
    a23 = np.einsum(
        "n,r,nra,ab,bsm,s,m->",
        e2,
        e3,
        v3(k2, k3, -p23),
        G,
        v3(p23, k4, k1),
        e4,
        e1,
    ) / mdot(p23, p23)

    # Quartic color-ordered tensor:
    # 2 g^{mu rho} g^{nu sigma} - g^{mu nu} g^{rho sigma}
    # - g^{mu sigma} g^{nu rho}.
    c12_34 = mdot(eps[0], eps[1]) * mdot(eps[2], eps[3])
    c13_24 = mdot(eps[0], eps[2]) * mdot(eps[1], eps[3])
    c14_23 = mdot(eps[0], eps[3]) * mdot(eps[1], eps[2])
    return a12, a23, c12_34, c13_24, c14_23


def amplitude(ks: list[np.ndarray], eps: list[np.ndarray]) -> float:
    a12, a23, c12_34, c13_24, c14_23 = four_gluon_terms(ks, eps)
    return float(a12 + a23 - c12_34 + 2.0 * c13_24 - c14_23)


def current_for_leg(ks: list[np.ndarray], eps: list[np.ndarray], leg: int) -> np.ndarray:
    """Recover the contravariant current J^A exposed at one external leg."""
    j = np.zeros(5)
    for a in range(5):
        basis = np.zeros(5)
        basis[a] = 1.0
        trial = list(eps)
        trial[leg] = basis
        # A = epsilon_A J^A = g_aa epsilon^a J^a for a basis vector.
        j[a] = G[a, a] * amplitude(ks, trial)
    return j


def kinematics(e: float, theta: float) -> list[np.ndarray]:
    """Two KK massive legs (1,2) and two genuine 4D gluon legs (3,4)."""
    k1 = np.array([-e, 0.0, 0.0, 0.0, -e])
    k2 = np.array([-e, 0.0, 0.0, 0.0, +e])
    k3 = np.array([+e, e * np.cos(theta), e * np.sin(theta), 0.0, 0.0])
    k4 = np.array([+e, -e * np.cos(theta), -e * np.sin(theta), 0.0, 0.0])
    return [k1, k2, k3, k4]


def main() -> None:
    rng = np.random.default_rng(20260826)
    ward_max = 0.0
    relation_max = 0.0
    longitudinal_values: list[float] = []

    for trial in range(20):
        e = 0.9 + 0.07 * trial
        theta = 0.23 + 0.041 * trial
        ks = kinematics(e, theta)

        assert np.linalg.norm(sum(ks)) < 1e-12
        assert max(abs(mdot(k, k)) for k in ks) < 1e-12

        eps = [transverse(k, rng) for k in ks]
        assert max(abs(mdot(k, ep)) for k, ep in zip(ks, eps)) < 1e-11

        # Full five-dimensional Ward identities.
        for leg in range(4):
            trial_eps = list(eps)
            trial_eps[leg] = ks[leg]
            ward_max = max(ward_max, abs(amplitude(ks, trial_eps)))

        # Expose the current of the first KK leg.  The other three legs remain
        # physical/transverse.
        j = current_for_leg(ks, eps, 0)
        k = ks[0]
        p4 = k[:4]
        j4 = j[:4]
        kappa = k[4]

        p_dot_j4 = float(p4 @ G[:4, :4] @ j4)
        relation = p_dot_j4 - kappa * j[4]
        relation_max = max(relation_max, abs(relation))
        longitudinal_values.append(abs(p_dot_j4))

        # Full 5D Ward identity is the same relation.
        assert abs(mdot(k, j)) < 2e-10

    print(f"max full 5D Ward residual: {ward_max:.3e}")
    print(f"max p.J4 - kappa*J5 residual: {relation_max:.3e}")
    print(f"min |p.J4| over trials: {min(longitudinal_values):.6e}")
    print(f"max |p.J4| over trials: {max(longitudinal_values):.6e}")

    # The important negative result: four-dimensional transversality is not true.
    if max(longitudinal_values) < 1e-6:
        raise AssertionError("unexpectedly found only 4D-transverse currents")
    if ward_max > 2e-10 or relation_max > 2e-10:
        raise AssertionError("Ward reconstruction audit failed")

    print("PASS: 5D Ward identity holds while generic 4D current transversality fails.")
    print("PASS: longitudinal massive-projector term reconstructs the fifth-current piece.")


if __name__ == "__main__":
    main()
