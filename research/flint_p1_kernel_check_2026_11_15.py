#!/usr/bin/env python
"""Flint eleventh run (2026-11-15): pair-internal check of Heart's P1.

Heart Round 5 (GPPDiscovery5 3ee15a8) certified, in arb ball
arithmetic on the SPECTRAL side, the finite odd block at L = 0.5:
    B(phi_j, phi_k) on E_12 = span{ sqrt(2/l) sin(w_j x), w_j = 4 pi j },
    l = 1/4,  >=  beta I,   beta = 1.66823e-8 = 0.9 * lambda_min,
    lambda_min(M_12) = 1.8535938681e-8 (float, from her ball midpoints).
Her own cross-check recomputed five entries kernel-side by float
quadrature.  THIS script is Flint's independent code path: the full
12x12 block assembled in the KERNEL domain from the in-tree
definitions (GPPVerify KreinHPsi.lean, read firsthand):
    K_S(d) = poleTerm(d) - prime_S(d) + archTerm(d),
    poleTerm(d) = 2 * bigPsiI2^2 cosh(d/2),  bigPsiI2 = 16 sinh^2(1/4),
    prime_S(d)  = sum_{n in {2,3,4,5,7,8,9,11}} Lambda(n)/sqrt(n)
                    * [triG(log n + d) + triG(-log n + d)],
    triG(y)     = 2/3 - y^2 + |y|^3/2 (|y|<=1); (2-|y|)^3/6 (1<=|y|<=2),
    archTerm(d) = (1/pi) int_0^inf bigPsi(t)^2 cos(t d)
                    (Re psi(1/4 + i t/2) - log pi) dt,
    bigPsi(t)   = (sin(t/2)/(t/2))^2,
    G(s, u)     = (1/2)[K_S(s - u) - K_S(s + u)]   on [0, l]^2,
    M_jk        = int int phi_j(s) phi_k(u) G(s, u) ds du.
(On |d| <= L = 0.5 the in-tree primeTerm and prime_S coincide: every
n outside the set has triG(+-log n + d) = 0 identically there.)

Method: float64 numpy.  archTerm is evaluated by period-wise
Gauss-Legendre (64 nodes per 2 pi period, t <= 4000; tail majorant
16 (log T + 2.83) / (3 T^3) ~ 9.3e-10) on a d-grid of 2001 points in
[0, 0.5] and PCHIP-interpolated; the evaluator is certified in-run
against the banked accurate value archTerm(0) = -1.54293308589872
(Heart, 96 GL nodes/period, cross-checked by quadosc; Turbo audit
6888858952-era: two independent float64 routes give -1.542933085912).
The plain-quad value -1.542933059454 is the known artifact and must
NOT appear.  The (s, u) integrals use composite Gauss-Legendre
(12 panels x 10 nodes per axis).  Prime and pole parts are exact
closed forms.  No zero data anywhere.

Checks against Heart's filed M (midpoints, round5 JSON):
entries (1,1), (1,2), (3,4), (6,7), (12,12); full spectrum and
lambda_min; beta = 0.9 lambda_min arithmetic; LDL^T pivots of
M_filed - beta I all positive (downstream arithmetic check).

Epistemic status: NUMERICAL (float64, independent code path).
This is a spot verification of a NUMERICAL certificate, not a
re-certification: the arb ball certificate is Heart's.
"""

import json
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import PchipInterpolator
from scipy.special import digamma

T_START = time.time()
OUT = {}

# ------------------------------------------------------------ archTerm
gx, gw = leggauss(64)
periods = int(np.ceil(4000 / (2 * np.pi)))
edges = 2 * np.pi * np.arange(periods + 1)
centers = 0.5 * (edges[:-1] + edges[1:])
halves = 0.5 * (edges[1:] - edges[:-1])
# nodes per period: t = center + half * gx
T_NODES = (centers[:, None] + halves[:, None] * gx[None, :]).ravel()
T_W = (halves[:, None] * gw[None, :]).ravel()
BIGPSI2 = (np.sin(T_NODES / 2) / (T_NODES / 2)) ** 4     # bigPsi(t)^2
WGT = (digamma(0.25 + 0.5j * T_NODES).real - np.log(np.pi))
BASE = T_W * BIGPSI2 * WGT / np.pi                      # arch(d) = sum BASE * cos(t d)

dgrid = np.linspace(0.0, 0.5, 2001)
arch_grid = np.array([np.sum(BASE * np.cos(T_NODES * d)) for d in dgrid])
arch_interp = PchipInterpolator(dgrid, arch_grid)
ARCH0_BANKED = -1.54293308589872
OUT['arch0_mine'] = float(arch_grid[0])
OUT['arch0_banked'] = ARCH0_BANKED
OUT['arch0_absdiff'] = float(abs(arch_grid[0] - ARCH0_BANKED))
print('archTerm(0): mine', arch_grid[0], ' banked', ARCH0_BANKED,
      ' |diff|', OUT['arch0_absdiff'])
tail_major = 16 * (np.log(4000) + 2.83) / (3 * 4000.0 ** 3)
OUT['arch_tail_majorant'] = float(tail_major)
print('arch tail majorant beyond T=4000:', tail_major)

# ------------------------------------------------------------ K_S pieces
BIGPSI_I2 = 16 * np.sinh(0.25) ** 2
LAM = {2: np.log(2), 3: np.log(3), 4: np.log(2), 5: np.log(5),
       7: np.log(7), 8: np.log(2), 9: np.log(3), 11: np.log(11)}


def triG(y):
    a = np.abs(y)
    return np.where(a <= 1, 2 / 3 - y ** 2 + a ** 3 / 2,
                    np.where(a <= 2, (2 - a) ** 3 / 6, 0.0))


def prime_S(d):
    tot = np.zeros_like(np.asarray(d, dtype=float))
    for n, lam in LAM.items():
        ln = np.log(n)
        tot = tot + lam / np.sqrt(n) * (triG(ln + d) + triG(-ln + d))
    return tot


def K_S(d):
    return (2 * BIGPSI_I2 ** 2 * np.cosh(d / 2)
            - prime_S(d) + arch_interp(np.abs(d)))


# ------------------------------------------------------------ block assembly
ell = 0.25
panels = 12
pn, pw = leggauss(10)
pedges = np.linspace(0, ell, panels + 1)
pc = 0.5 * (pedges[:-1] + pedges[1:])
ph = 0.5 * (pedges[1:] - pedges[:-1])
S_NODES = (pc[:, None] + ph[:, None] * pn[None, :]).ravel()
S_W = (ph[:, None] * pw[None, :]).ravel()
M_NODES = len(S_NODES)
print('quadrature nodes per axis:', M_NODES)

Dminus = S_NODES[:, None] - S_NODES[None, :]
Dplus = S_NODES[:, None] + S_NODES[None, :]
Gmat = 0.5 * (K_S(Dminus) - K_S(Dplus))
WG = S_W[:, None] * Gmat * S_W[None, :]

PHI = np.array([np.sqrt(2 / ell) * np.sin(4 * np.pi * j * S_NODES)
                for j in range(1, 13)])                  # (12, M)
M_mine = PHI @ WG @ PHI.T
M_mine = 0.5 * (M_mine + M_mine.T)
ev_mine = np.linalg.eigvalsh(M_mine)
OUT['lambda_min_mine'] = float(ev_mine[0])
OUT['eigs_mine'] = [float(x) for x in ev_mine]
OUT['M_mine'] = [[float(x) for x in row] for row in M_mine]
print('lambda_min (mine, kernel domain):', ev_mine[0])
print('spectrum (mine):', ['%.4e' % x for x in ev_mine])

# ------------------------------------------------------------ filed values
filed = json.load(open('/tmp/heart/r5_N12.json'))
M_filed = np.array([[float(x) for x in row] for row in filed['M']])
ev_filed = np.linalg.eigvalsh(0.5 * (M_filed + M_filed.T))
OUT['lambda_min_filed_midpoints'] = float(ev_filed[0])
print('lambda_min (her filed midpoints):', ev_filed[0],
      ' her float value:', filed['lambda_min_float'])
spot = {}
for (j, k) in [(0, 0), (0, 1), (2, 3), (5, 6), (11, 11)]:
    spot[f'M{j+1}{k+1}'] = {
        'mine': float(M_mine[j, k]), 'filed': float(M_filed[j, k]),
        'absdiff': float(abs(M_mine[j, k] - M_filed[j, k])),
        'reldiff': float(abs(M_mine[j, k] - M_filed[j, k]) / abs(M_filed[j, k]))}
    print(f'M[{j+1},{k+1}] mine {M_mine[j,k]:.12e} filed {M_filed[j,k]:.12e} '
          f'rel {spot[f"M{j+1}{k+1}"]["reldiff"]:.3e}')
OUT['spot_entries'] = spot
OUT['lambda_min_absdiff'] = float(abs(ev_mine[0] - ev_filed[0]))

beta_exact = 0.9 * float(filed['lambda_min_float'])
OUT['beta_from_filed_lambda_min'] = beta_exact
OUT['beta_filed'] = float(filed['beta_certified'].strip('[]').split(' ')[0])
print('beta: 0.9 * her lambda_min =', beta_exact, ' filed beta =',
      OUT['beta_filed'], '(truncated down:', beta_exact > OUT['beta_filed'], ')')

# LDL^T pivots of M_filed - beta I (downstream arithmetic check)
A = M_filed - OUT['beta_filed'] * np.eye(12)
L = np.zeros((12, 12))
D = np.zeros(12)
for j in range(12):
    D[j] = A[j, j] - sum(L[j, m] ** 2 * D[m] for m in range(j))
    for i in range(j + 1, 12):
        L[i, j] = (A[i, j] - sum(L[i, m] * L[j, m] * D[m] for m in range(j))) / D[j]
OUT['ldl_pivots_Mfiled_minus_betaI'] = [float(x) for x in D]
print('LDL pivots of M_filed - beta I all positive:',
      bool(np.all(D > 0)), ' min pivot', float(np.min(D)))

OUT['runtime_seconds'] = time.time() - T_START
with open('/home/hatch/workspace/gppdiscovery3/research/'
          'flint_p1_kernel_check_2026_11_15_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('done', OUT['runtime_seconds'])
