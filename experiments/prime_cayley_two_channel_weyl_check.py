"""Finite numerical audit of the prime-energy Cayley / two-channel Weyl construction.

Discovery-only: verifies finite truncations of exact identities recorded in
GPP-bridge research/codex/2026-10-02_prime_energy_cayley_two_channel_weyl.md.
No RH claim.
"""

import numpy as np
from math import sqrt

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]


def S_from_p(p: int) -> np.ndarray:
    q = p ** -0.5
    b = sqrt(1.0 - 1.0 / p)
    return np.array([[q, b], [b, -q]], dtype=float)


def S_from_mu(p: int) -> np.ndarray:
    mu = 2.0 / sqrt(p - 1.0)
    return np.array([[mu, 2.0], [2.0, -mu]]) / sqrt(mu * mu + 4.0)


TAU = np.array([[0.0, -1.0], [1.0, 0.0]])


def weyl_matrix(z: complex, primes=PRIMES) -> np.ndarray:
    out = np.zeros((2, 2), dtype=complex)
    for p in primes:
        v = np.array([p ** -0.25, p ** -0.5 / sqrt(2.0)], dtype=complex)
        vv = np.outer(v, np.conjugate(v))
        out += vv * (1.0 / (p - z) - p / (1.0 + p * p))
    return out


def schur_cayley(M: np.ndarray) -> np.ndarray:
    I = np.eye(M.shape[0], dtype=complex)
    return (M - 1j * I) @ np.linalg.inv(M + 1j * I)


def main():
    tol = 2e-12

    for p in PRIMES:
        S = S_from_p(p)
        Sm = S_from_mu(p)
        assert np.linalg.norm(S - Sm) < tol
        assert np.linalg.norm(S @ S - np.eye(2)) < tol
        assert np.linalg.norm(S @ TAU + TAU @ S) < tol

    phi = (1.0 + sqrt(5.0)) / 2.0
    c5 = (sqrt(5.0) - 1.0) / (sqrt(5.0) + 1.0)
    assert abs(c5 - phi ** -2) < tol

    test_points = [
        -3.0 + 0.2j,
        0.0 + 0.5j,
        4.5 + 0.7j,
        20.0 + 3.0j,
    ]
    for z in test_points:
        M = weyl_matrix(z)
        ImM = (M - M.conj().T) / (2j)
        eig_im = np.linalg.eigvalsh(ImM)
        assert eig_im.min() > -tol

        Sigma = schur_cayley(M)
        defect = np.eye(2) - Sigma.conj().T @ Sigma
        eig_defect = np.linalg.eigvalsh((defect + defect.conj().T) / 2)
        assert eig_defect.min() > -5e-11
        assert np.linalg.svd(Sigma, compute_uv=False).max() <= 1.0 + 5e-11

    Mi = weyl_matrix(1j)
    assert np.linalg.norm(Mi.real) < 5e-12

    print("finite prime Cayley/Weyl audit passed")
    print("p=5 c =", c5, "phi^-2 =", phi ** -2)
    for z in test_points:
        M = weyl_matrix(z)
        Sigma = schur_cayley(M)
        print(
            "z=", z,
            "min eig Im M=", np.linalg.eigvalsh((M-M.conj().T)/(2j)).min(),
            "||Sigma||=", np.linalg.svd(Sigma, compute_uv=False).max(),
        )


if __name__ == "__main__":
    main()
