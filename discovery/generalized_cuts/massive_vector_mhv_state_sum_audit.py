#!/usr/bin/env python3
"""Audit the adjacent-MHV massive-vector state sum in the 5D embedding.

Codex/GPT Golden Physics discovery track, 2026-08-26.

We evaluate the COMPLETE color-ordered Yang-Mills four-point tree (both cubic
exchange channels plus the quartic contact term) for

    V_m(K1) V_m(K2) g^h(k3) g^h(k4),   h = +/-

using the 5D massless embedding K=(p,kappa), with 4D external gluons.  In the
centre-of-mass frame the three physical massive-vector states are the three
ordinary spatial unit vectors.  Generic 2->2 kinematics are parameterized by
energy E and scattering angle theta.

In the normalization of the standard color-ordered vertices used here, the
same-helicity tree matrix is found to be

    A_h^{ab} = - delta^{ab},  a,b=1,2,3,

for both h=+ and h=-, independently of E and theta.  Embedding one dimension
higher and choosing the sixth polarization gives the real-adjoint-scalar tree

    A_h^(S) = +1.

Consequently, contracting the two equal-helicity tree factors across a
(-- | ++) two-particle cut gives the state-count ratio

    C^(V_m) / C^(S) = Tr(I_3) = 3,

and dimensional reconstruction then predicts in this sector

    C^(4) = C^(V_m) - C^(S) = 2 C^(S),
    C^(D_s) = (D_s - 2) C^(S).

This script audits the tree/state algebra only.  Overall coupling, color,
spinor-helicity prefactors, cut orientation, and loop-measure normalization are
not fixed here; the ratio is the invariant target to map into the paper's cut
normalization before formal promotion.
"""

from __future__ import annotations

import numpy as np


def metric(d: int) -> np.ndarray:
    return np.diag([1.0] + [-1.0] * (d - 1)).astype(complex)


def mdot(g: np.ndarray, a: np.ndarray, b: np.ndarray) -> complex:
    return a @ g @ b


def v3(g: np.ndarray, p: np.ndarray, q: np.ndarray, r: np.ndarray) -> np.ndarray:
    d = len(p)
    out = np.zeros((d, d, d), dtype=complex)
    for mu in range(d):
        for nu in range(d):
            for rho in range(d):
                out[mu, nu, rho] = (
                    g[mu, nu] * (p - q)[rho]
                    + g[nu, rho] * (q - r)[mu]
                    + g[rho, mu] * (r - p)[nu]
                )
    return out


def amplitude(ks: list[np.ndarray], eps: list[np.ndarray]) -> complex:
    """Complete color-ordered 4-gluon tree in mostly-minus signature."""
    d = len(ks[0])
    g = metric(d)
    k1, k2, k3, k4 = ks
    e1, e2, e3, e4 = [g @ e for e in eps]

    p12 = k1 + k2
    a12 = np.einsum(
        "m,n,mna,ab,brs,r,s->",
        e1, e2, v3(g, k1, k2, -p12), g,
        v3(g, p12, k3, k4), e3, e4,
    ) / mdot(g, p12, p12)

    p23 = k2 + k3
    a23 = np.einsum(
        "n,r,nra,ab,bsm,s,m->",
        e2, e3, v3(g, k2, k3, -p23), g,
        v3(g, p23, k4, k1), e4, e1,
    ) / mdot(g, p23, p23)

    c12_34 = mdot(g, eps[0], eps[1]) * mdot(g, eps[2], eps[3])
    c13_24 = mdot(g, eps[0], eps[2]) * mdot(g, eps[1], eps[3])
    c14_23 = mdot(g, eps[0], eps[3]) * mdot(g, eps[1], eps[2])

    # Color-ordered quartic tensor:
    # 2 g^{mu rho} g^{nu sigma} - g^{mu nu} g^{rho sigma}
    # - g^{mu sigma} g^{nu rho}.
    return a12 + a23 - c12_34 + 2.0 * c13_24 - c14_23


def kinematics(e: float, theta: float, d: int = 5) -> list[np.ndarray]:
    """All-outgoing 5D null kinematics with two opposite KK charges."""
    base = [
        np.array([-e, 0.0, 0.0, 0.0, -e], dtype=complex),
        np.array([-e, 0.0, 0.0, 0.0, +e], dtype=complex),
        np.array([+e, e * np.cos(theta), e * np.sin(theta), 0.0, 0.0], dtype=complex),
        np.array([+e, -e * np.cos(theta), -e * np.sin(theta), 0.0, 0.0], dtype=complex),
    ]
    if d == 5:
        return base
    if d == 6:
        return [np.concatenate([k, np.zeros(1, dtype=complex)]) for k in base]
    raise ValueError("only d=5 or d=6 used in this audit")


def gluon_helicity(theta: float, leg: int, h: int, d: int = 5) -> np.ndarray:
    """A transverse circular-polarization basis for k3 or k4."""
    if leg == 3:
        e_plane = np.array([0.0, -np.sin(theta), np.cos(theta), 0.0, 0.0], dtype=complex)
    elif leg == 4:
        e_plane = np.array([0.0, +np.sin(theta), -np.cos(theta), 0.0, 0.0], dtype=complex)
    else:
        raise ValueError("leg must be 3 or 4")
    e_z = np.array([0.0, 0.0, 0.0, 1.0, 0.0], dtype=complex)
    eps = (e_plane + 1j * h * e_z) / np.sqrt(2.0)
    if d == 6:
        eps = np.concatenate([eps, np.zeros(1, dtype=complex)])
    return eps


def massive_basis() -> list[np.ndarray]:
    return [
        np.array([0, 1, 0, 0, 0], dtype=complex),
        np.array([0, 0, 1, 0, 0], dtype=complex),
        np.array([0, 0, 0, 1, 0], dtype=complex),
    ]


def vector_tree_matrix(e: float, theta: float, h: int) -> np.ndarray:
    ks = kinematics(e, theta, 5)
    ep3 = gluon_helicity(theta, 3, h, 5)
    ep4 = gluon_helicity(theta, 4, h, 5)
    basis = massive_basis()
    return np.array(
        [[amplitude(ks, [ea, eb, ep3, ep4]) for eb in basis] for ea in basis],
        dtype=complex,
    )


def scalar_tree(e: float, theta: float, h: int) -> complex:
    ks = kinematics(e, theta, 6)
    ep3 = gluon_helicity(theta, 3, h, 6)
    ep4 = gluon_helicity(theta, 4, h, 6)
    scalar = np.array([0, 0, 0, 0, 0, 1], dtype=complex)
    return amplitude(ks, [scalar, scalar, ep3, ep4])


def main() -> None:
    max_matrix_residual = 0.0
    max_scalar_residual = 0.0
    max_ratio_residual = 0.0

    # Avoid forward/backward singular configurations; otherwise scan generic 2->2 data.
    energies = [0.73, 1.0, 1.37, 2.11, 4.2]
    angles = [0.19, 0.41, 0.77, 1.13, 1.49, 2.03, 2.61]

    for e in energies:
        for theta in angles:
            for h in (-1, +1):
                M = vector_tree_matrix(e, theta, h)
                S = scalar_tree(e, theta, h)

                matrix_residual = np.max(np.abs(M + np.eye(3)))
                scalar_residual = abs(S - 1.0)
                max_matrix_residual = max(max_matrix_residual, float(matrix_residual))
                max_scalar_residual = max(max_scalar_residual, float(scalar_residual))

            Mminus = vector_tree_matrix(e, theta, -1)
            Mplus = vector_tree_matrix(e, theta, +1)
            Sminus = scalar_tree(e, theta, -1)
            Splus = scalar_tree(e, theta, +1)

            # Canonical internal-state contraction.  With M_-=M_+=-I this is 3.
            cv = np.einsum("ab,ab->", Mminus, Mplus)
            cs = Sminus * Splus
            max_ratio_residual = max(max_ratio_residual, float(abs(cv - 3.0 * cs)))

    print(f"max ||A_h + I_3||_inf: {max_matrix_residual:.3e}")
    print(f"max |A_h^(S) - 1|:     {max_scalar_residual:.3e}")
    print(f"max |C_V - 3 C_S|:     {max_ratio_residual:.3e}")

    tol = 2e-11
    if max_matrix_residual > tol:
        raise AssertionError("massive-vector same-helicity tree is not -I_3")
    if max_scalar_residual > tol:
        raise AssertionError("extra-dimensional scalar same-helicity tree is not +1")
    if max_ratio_residual > tol:
        raise AssertionError("massive-vector/scalar sewn state-count ratio is not 3")

    print("PASS: complete same-helicity vector tree matrix = -I_3.")
    print("PASS: complete extra-dimensional scalar tree = +1.")
    print("PASS: adjacent (--|++) state sewing gives C^(V_m)=3 C^(S).")
    print("TARGET: after normalization mapping, C^(D_s)=(D_s-2) C^(S) in this sector.")


if __name__ == "__main__":
    main()
