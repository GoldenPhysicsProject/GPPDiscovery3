#!/usr/bin/env python
"""Flint tenth run (2026-11-14): candidate (b) computations at the
clustered L = 0.5, n = 6 configuration (x_i in [-0.025, 0.025]).

PART 1 -- R2 leg (i), spectral-face form.  The kernel-face projected
lambda_min (anchor +8.88e-7, kill arm / F_theta control) is recomputed
at dps 50, then reproduced from the SPECTRAL (state) face alone via
identity C1: the matrix
    M_spec[i,j] = sum_{all rho} K(rho) K(1-rho) e^{-(rho-1/2)(x_i-x_j)},
    K(rho) = (sinh((rho-1/2)/2)/((rho-1/2)/2))^2
            = (sin(gamma/2)/(gamma/2))^2  on the line,
plus the pole boundary matrix P of identity (I3)
(c^T P c = Ftilde_c(0)^2 + Ftilde_c(1)^2).  No Kpsi kernel value is used
on the state face.  Platt zeros |gamma| <= 3000, analytic tail bound.

PART 2 -- Heart's C1 (her Round 1, GPPDiscovery5 41dcd57; her conjecture,
registered there before this run): the parity split.  The configuration
is symmetric under reversal (x -> -x), so D block-diagonalises:
  even block: vectors even under reversal, intersected with ker l_c
              (the pole term vanishes identically: sum c_i sinh(x_i/2)=0);
  odd block:  vectors odd under reversal (ker l_c automatic:
              sum c_i cosh(x_i/2) = 0), full D = arch - prime - pole.
Her prediction: the odd face carries the small margin (+8.88e-7).
Both faces computed on the kernel face (dps 50) and on the state face
(M_spec + P) as a cross-check.

Kernel-face conventions copied from flint_ftheta_control_2026_11_13.py
(arch_mp: (sin(t/2)/(t/2))^4 digamma kernel; prime_mp: triG prime sums,
theta = False).  Epistemic status: NUMERICAL.
"""

import json
import sys
import time

sys.path.insert(0, '/home/hatch/workspace/zeros-platt')
sys.path.insert(0, '/home/hatch/workspace/gppdiscovery3')

import numpy as np
import mpmath
from mpmath import mp, mpf
mp.dps = 50

T_START = time.time()
OUT = {}

# ------------------------------------------------------------ configuration
import importlib.util
spec = importlib.util.spec_from_file_location(
    'scout', '/home/hatch/workspace/gppdiscovery3/research/'
             'flint_archdom_scout_2026_11_12.py')
scout = importlib.util.module_from_spec(spec)
import types
# the scout runs a main() only under __main__ guard check: load safely
src = open(scout.__dict__.get('__file__', '') or
           '/home/hatch/workspace/gppdiscovery3/research/'
           'flint_archdom_scout_2026_11_12.py').read()
NS = {}
exec(compile(src.split("if __name__")[0], 'scout-defs', 'exec'), NS)
configs = NS['configs']
x = np.array(configs(0.5, 6)['clustered0.1L'], dtype=float)
n = len(x)
print('config x =', x)
OUT['x'] = [str(v) for v in x]

# ------------------------------------------------------------ kernel face (mp)
def triG_mp(y):
    ay = abs(y)
    if ay <= 1:
        return mpf(2) / 3 - y * y + ay ** 3 / 2
    if ay <= 2:
        return (2 - ay) ** 3 / 6
    return mpf(0)


def prime_powers(nmax):
    out = []
    sieve = [True] * (nmax + 1)
    for p in range(2, nmax + 1):
        if sieve[p]:
            for j in range(p * p, nmax + 1, p):
                sieve[j] = False
            pk = p
            while pk <= nmax:
                out.append((pk, mp.log(p)))
                if pk > nmax // p:
                    break
                pk *= p
    return sorted(out)


PP = prime_powers(64)


def prime_mp(d):
    d = mpf(d)
    lim = mp.e ** (2 + abs(d))
    s = mpf(0)
    for nn, lp in PP:
        if nn > lim:
            break
        ln = mp.log(nn)
        s += lp / mp.sqrt(nn) * (triG_mp(ln + d) + triG_mp(-ln + d))
    return s


def arch_mp(d):
    d = mpf(d)
    f = lambda t: (mpf(1) if t == 0 else (mp.sin(t / 2) / (t / 2)) ** 2) ** 2 \
        * mp.cos(t * d) * (mp.re(mp.digamma(mp.mpc(mpf('0.25'), t / 2))) - mp.log(mp.pi))
    return (1 / mp.pi) * mp.quad(f, [0, mp.inf])


bI2 = 16 * mp.sinh(mpf('0.25')) ** 2
Dmp = [[mpf(0)] * n for _ in range(n)]
Amp_cache, Pth_cache = {}, {}
for i in range(n):
    for j in range(n):
        d = abs(mpf(str(x[i])) - mpf(str(x[j])))
        key = mp.nstr(d, 12)
        if key not in Amp_cache:
            Amp_cache[key] = arch_mp(d)
            Pth_cache[key] = prime_mp(d)
        Dmp[i][j] = (Amp_cache[key] - Pth_cache[key]
                     - 2 * bI2 ** 2 * mp.sinh(mpf(str(x[i])) / 2)
                     * mp.sinh(mpf(str(x[j])) / 2))
D = np.array([[float(Dmp[i][j]) for j in range(n)] for i in range(n)])
print('kernel-face D built; D[0,0] =', mp.nstr(Dmp[0][0], 15))
OUT['D_kernel'] = [[mp.nstr(Dmp[i][j], 15) for j in range(n)] for i in range(n)]


def proj_lammin(M, x):
    ell = np.cosh(x / 2)
    _, _, Vt = np.linalg.svd(ell.reshape(1, -1))
    Z = Vt[1:].T
    Mp = Z.T @ M @ Z
    w, v = np.linalg.eigh(Mp)
    return w[0], Z @ v[:, 0], Z


lam_D, cstar, Zbasis = proj_lammin(D, x)
print('kernel face: projected lambda_min =', lam_D)
OUT['lambda_min_kernel'] = lam_D
# certified quadratic form at cstar, dps 50
cs = [mpf(str(v)) for v in cstar]
qD = sum(cs[i] * Dmp[i][j] * cs[j] for i in range(n) for j in range(n))
cnorm = mp.sqrt(sum(v * v for v in cs))
print('Q_D(c*) dps50 =', mp.nstr(qD, 15), ' |c*| =', mp.nstr(cnorm, 10))
OUT['Q_D_cstar'] = mp.nstr(qD, 18)
OUT['cstar'] = [mp.nstr(v, 15) for v in cs]

# ------------------------------------------------------------ state face
from platt_sequential import zeros_in_range
T = 3000
M = np.zeros((n, n))
for g in zeros_in_range(0, T):
    gf = float(g)
    kern = (np.sin(gf / 2) / (gf / 2)) ** 4
    for i in range(n):
        for j in range(n):
            M[i, j] += 2 * kern * np.cos(gf * (x[i] - x[j]))
tail = (16 / np.pi) * (np.log(T) + 1 / 3) / T ** 3
print('state face M built; entry tail bound ~', tail)
OUT['M_tail_bound'] = tail
ep = np.exp(x / 2); em = np.exp(-x / 2)
P = float(bI2) ** 2 * (np.outer(ep, ep) + np.outer(em, em))
for name, Mx in (('M_spec', M), ('M_spec+P', M + P)):
    lam, _, _ = proj_lammin(Mx, x)
    q = float(cstar @ Mx @ cstar)
    print(f'{name}: projected lambda_min = {lam:.6e}   Q(c*) = {q:.6e}')
    OUT['lambda_min_' + name] = lam
    OUT['Q_' + name + '_cstar'] = q

# ------------------------------------------------------------ parity split (C1)
rev = np.array([[1.0 if j == n - 1 - i else 0.0 for j in range(n)]
                for i in range(n)])
# even/odd bases under reversal
Ebasis, Obasis = [], []
for i in range(n // 2):
    e = np.zeros(n); e[i] = 1 / np.sqrt(2); e[n - 1 - i] = 1 / np.sqrt(2)
    o = np.zeros(n); o[i] = 1 / np.sqrt(2); o[n - 1 - i] = -1 / np.sqrt(2)
    Ebasis.append(e); Obasis.append(o)
E = np.array(Ebasis).T; O = np.array(Obasis).T
# sanity: D commutes with reversal; pole vanishes on even; constraint auto on odd
print('|D rev - rev D|_max =', np.abs(D @ rev - rev @ D).max())
ellv = np.cosh(x / 2); sinhv = np.sinh(x / 2)
print('max |ell . odd basis| =', np.abs(ellv @ O).max(),
      ' max |sinh . even basis| =', np.abs(sinhv @ E).max())
OUT['parity_sanity'] = {'commutator': float(np.abs(D @ rev - rev @ D).max()),
                        'ell_on_odd': float(np.abs(ellv @ O).max()),
                        'sinh_on_even': float(np.abs(sinhv @ E).max())}


def block_lammin(Mx):
    # odd block: unconstrained
    Mo = O.T @ Mx @ O
    lam_o = np.linalg.eigvalsh(Mo)[0]
    # even block: intersect ker l_c within even coords
    Me = E.T @ Mx @ E
    elle = E.T @ ellv
    _, _, Vt = np.linalg.svd(elle.reshape(1, -1))
    Ze = Vt[1:].T
    lam_e = np.linalg.eigvalsh(Ze.T @ Me @ Ze)[0]
    return lam_e, lam_o


for name, Mx in (('kernel', D), ('state', M + P)):
    lam_e, lam_o = block_lammin(Mx)
    print(f'C1 {name}: lambda_min even block = {lam_e:.6e}   '
          f'odd block = {lam_o:.6e}')
    OUT[f'C1_{name}'] = {'lambda_min_even': lam_e, 'lambda_min_odd': lam_o}

OUT['runtime_seconds'] = time.time() - T_START
with open('/home/hatch/workspace/gppdiscovery3/research/'
          'flint_r2_leg1_spectral_2026_11_14_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('done', OUT['runtime_seconds'])
