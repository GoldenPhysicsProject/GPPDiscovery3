#!/usr/bin/env python3
"""Audit the 2-3-5-7 arithmetic-family incidence hypothesis.

Exploratory calculation: exact graph identities are separated from
phenomenological tests. No Standard-Model derivation is claimed.

Core objects:
  B    = incidence of 3--5--7
  L_nu = B diag(a,b) B^T
  L_ch = principal Dirichlet minor of the 2--3--5--7 Laplacian
"""

from __future__ import annotations

import math
import numpy as np


def neutral_laplacian(a: float, b: float) -> np.ndarray:
    return np.array([
        [a, -a, 0.0],
        [-a, a + b, -b],
        [0.0, -b, b],
    ], dtype=float)


def charged_dirichlet(c: float, a: float, b: float) -> np.ndarray:
    # Pin the boundary node p=2 in the chain 2--3--5--7.
    return np.array([
        [c + a, -a, 0.0],
        [-a, a + b, -b],
        [0.0, -b, b],
    ], dtype=float)


def edge_ratio_for_target(q: float) -> tuple[float, float]:
    """Solve lambda_-/lambda_+ = q for a=1, b=t."""
    def f(t: float) -> float:
        d = math.sqrt(1.0 - t + t * t)
        return ((1.0 + t) - d) / ((1.0 + t) + d) - q

    def bisect(lo: float, hi: float) -> float:
        flo = f(lo)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            fm = f(mid)
            if abs(fm) < 1e-15:
                return mid
            if flo * fm <= 0:
                hi = mid
            else:
                lo, flo = mid, fm
        return 0.5 * (lo + hi)

    t1 = bisect(1e-8, 1.0)
    return t1, 1.0 / t1


def mixing_angles_from_abs(V: np.ndarray) -> tuple[float, float, float]:
    s13 = min(1.0, abs(V[0, 2]))
    c13 = math.sqrt(max(0.0, 1.0 - s13 * s13))
    s12 = min(1.0, abs(V[0, 1]) / c13)
    s23 = min(1.0, abs(V[1, 2]) / c13)
    return tuple(math.degrees(math.asin(x)) for x in (s12, s23, s13))


def main() -> None:
    B = np.array([[1.0, 0.0], [-1.0, 1.0], [0.0, -1.0]])
    print("B^T B =")
    print(B.T @ B)
    print("ker(B^T) representative:", np.array([1.0, 1.0, 1.0]))
    print()

    # Approximate 2025 PDG oscillation scales; only their ratio is used.
    dm21 = 7.49e-5
    dm3l = 2.50e-3
    target_mratio = math.sqrt(dm21 / dm3l)
    print("one-massless normal-order target m_light/m_heavy ~", target_mratio)
    print("required b/a branches for pure two-edge Laplacian:",
          edge_ratio_for_target(target_mratio))
    print()

    w = {p: math.log(p) / math.sqrt(p) for p in (2, 3, 5, 7)}
    q = {p: 1.0 / math.sqrt(p) for p in (2, 3, 5, 7)}
    kappa = {p: math.atanh(q[p]) for p in (2, 3, 5, 7)}
    h35 = 0.5 * math.log(5.0 / 3.0)
    h57 = 0.5 * math.log(7.0 / 5.0)

    schemes = {
        "equal": (1.0, 1.0),
        "Euler geometric": (math.sqrt(w[3] * w[5]), math.sqrt(w[5] * w[7])),
        "TFD-q geometric": (math.sqrt(q[3] * q[5]), math.sqrt(q[5] * q[7])),
        "kappa geometric": (math.sqrt(kappa[3] * kappa[5]),
                            math.sqrt(kappa[5] * kappa[7])),
        "log half-width": (h35, h57),
        "inverse log half-width": (1.0 / h35, 1.0 / h57),
    }

    print("neutral bare-path audit")
    for name, (a, b) in schemes.items():
        vals = np.linalg.eigvalsh(neutral_laplacian(a, b))
        r = vals[1] / vals[2]
        print(f"{name:24s} b/a={b/a: .9f} "
              f"lambda-/lambda+={r: .9f} "
              f"squared-ratio={r*r: .9f}")
    print()

    # Concrete no-free-parameter proxy: geometric Euler edge conductances,
    # including the exceptional p=2--3 boundary edge.
    a = math.sqrt(w[3] * w[5])
    b = math.sqrt(w[5] * w[7])
    c = math.sqrt(w[2] * w[3])
    Ln = neutral_laplacian(a, b)
    Lc = charged_dirichlet(c, a, b)
    en, Un = np.linalg.eigh(Ln)
    ec, Uc = np.linalg.eigh(Lc)
    Vproxy = Uc.T @ Un
    print("Euler-geometric neutral eigenvalues:", en)
    print("Euler-geometric charged Dirichlet eigenvalues:", ec)
    print("det(L_ch) numeric:", np.linalg.det(Lc), "abc:", a * b * c)
    print("|U_ch^T U_nu| proxy:")
    print(np.abs(Vproxy))
    print("proxy angles theta12,theta23,theta13 [deg]:",
          mixing_angles_from_abs(Vproxy))
    print()

    print("neutral zero mode (normalized, up to sign):",
          Un[:, np.argmin(np.abs(en))])
    print("A direct prime-family = (e,mu,tau) flavor identification would")
    print("therefore make the massless column democratic. Current data do not")
    print("have a democratic m1 (normal) or m3 (inverted) column; a nontrivial")
    print("charged-lepton rotation or weighted incidence map is required.")
    print()

    # General edge covariance C=[[a,x],[x,b]]. The two nonzero eigenvalues
    # have trace 2(a+b-x) and product 3(ab-x^2).
    def ratio_with_cross(x: float) -> float:
        T = 2.0 * (a + b - x)
        P = 3.0 * (a * b - x * x)
        disc = T * T - 4.0 * P
        return (T - math.sqrt(max(0.0, disc))) / (T + math.sqrt(max(0.0, disc)))

    sb = math.sqrt(a * b)
    roots = []
    grid = np.linspace(-0.999999 * sb, 0.999999 * sb, 20001)
    x0 = grid[0]
    y0 = ratio_with_cross(x0) - target_mratio
    for x1 in grid[1:]:
        y1 = ratio_with_cross(x1) - target_mratio
        if y0 * y1 < 0:
            lo, hi = x0, x1
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                ym = ratio_with_cross(mid) - target_mratio
                ylo = ratio_with_cross(lo) - target_mratio
                if ylo * ym <= 0:
                    hi = mid
                else:
                    lo = mid
            roots.append(0.5 * (lo + hi))
        x0, y0 = x1, y1

    print("edge-space cross covariance needed for observed one-massless ratio")
    for x in roots:
        print(f"x={x:.12f}, x/sqrt(ab)={x/sb:.12f}, "
              f"ratio={ratio_with_cross(x):.12f}")
    print()

    # Charged-lepton hierarchy falsifier for the same bare Dirichlet matrix.
    me = 0.51099895
    mmu = 105.6583755
    mtau = 1776.86
    observed = np.array([me, mmu, mtau]) / mtau
    proxy = ec / ec[-1]
    print("charged-lepton normalized masses [e,mu,tau]:", observed)
    print("bare Dirichlet normalized eigenvalues:", proxy)
    print("The boundary-lift mechanism fixes rank only; it does not generate")
    print("the charged-fermion hierarchy.")
    
    
    # Charge-selected p=2 boundary diagnostic.
    # Replacing c by c q_electric^2 gives det = a*b*c*q^2 and is even under q -> -q.
    charges = {"nu": 0.0, "e": -1.0, "u": 2.0 / 3.0, "d": -1.0 / 3.0}
    print()
    print("charge-selected p=2 boundary (same bare a,b,c; rank diagnostic only)")
    for name, qel in charges.items():
        Mq = charged_dirichlet(c * qel * qel, a, b)
        vals = np.linalg.eigvalsh(Mq)
        print(f"{name:2s} q={qel: .6f} det={np.linalg.det(Mq): .12e} eig={vals}")
    print("This selector explains neutral-vs-charged rank only; it is not a mass-hierarchy fit.")
    


if __name__ == "__main__":
    main()
