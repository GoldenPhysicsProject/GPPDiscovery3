#!/usr/bin/env python
"""Flint twelfth run (2026-11-15): the omega^-3 tail bound (b) for the
OddFace(0.5) closure attempt -- arm two, in parallel with Heart's (a).

CONTEXT.  Heart's Round-4 skeleton (GPPDiscovery5 0fbae22) splits
OddFace(0.5) into P1 (finite block E_12, closed Round 5, beta =
1.66823e-8, NUMERICAL interval-certified) + P2/P3; her Round 6
(2d3506b) corrected P2/P3: the high block is dominated by the boundary
rank-one term g*w_j*w_k (g = G(l,l) = 1.23787e-4, l = 1/4), scalar P3
is KILLED, the correct object is the bordered Schur complement
M11 - uu^T/g >= beta0 I with beta0 = 1.30163e-8 (NUMERICAL,
interval-certified), and the named GAP is: (a) arb-certify the
bordered Schur on Y = E_{12+K}, K ~ 60, at beta1; (b) bound the
omega^-3 tail loss from modes > 12+K below beta1.  THIS script is (b),
from Flint's own machinery.

OBJECT (exact reduction, Heart Round 5, CHECKED):
    B(phi_j, phi_k) = (2/l){ (1/pi) int_0^inf Psi(t)^2 V(t) S_j S_k dt
                             - 2*bigPsiI2^2 sigma_j sigma_k },
    l = 1/4, w_j = 4 pi j, Psi(t)^2 = (sin(t/2)/(t/2))^4,
    V(t) = Re psi(1/4 + i t/2) - log pi - 2 P(t),
    P(t) = sum_{n in {2,3,4,5,7,8,9,11}} Lambda(n) n^{-1/2} cos(t log n),
    S_j(t) = w_j l sinc((t-w_j) l)/(t + w_j)   [sinc(x) = sin x / x],
    sigma_j = -(-1)^j w_j sinh(l/2)/(1/4 + w_j^2),
    bigPsiI2 = 16 sinh^2(1/4).
Kernel side (in-tree, KreinHPsi.lean, as in run 11):
    K_S(d) = 2*bigPsiI2^2 cosh(d/2) - prime_S(d) + archTerm(d),
    G(s,u) = (1/2)[K_S(s-u) - K_S(s+u)],  b2(u) := d^2_s G(l,u)
           = (1/2)[K_S''(l-u) - K_S''(l+u)].

WHAT IS COMPUTED.
  1. Kernel constants: g = G(l,l) = (1/2)[K_S(0) - K_S(2l)] (kernel
     domain); u_j = B(phi_j, delta_l) for j <= 812 via the exact
     integration-by-parts identity u_j = g w_j - sqrt(2/l) w_j^{-2}
     int_0^l sin(w_j u) b2(u) du; d_i = <phi_i, b2> (i <= 12);
     h(l) := b2(l) = (1/2)[K_S''(0) - K_S''(2l)]; ||b2||_1, ||b2||_2,
     sup|b2|; C3 = int_0^l ||d^3_s G(s,.)||_2 ds; the R6.1 residual
     constant C_R = sqrt(2/l)(||b2||_2 + C3), so that
     ||(G phi_j) - w_j G(.,l)||_2 <= C_R w_j^{-3}.
  2. The full Gram M on modes 1..812 by period-wise Gauss-Legendre
     (48 nodes per 2 pi period, T = w_812 + 400; tail majorant
     ~5e-17 per entry).  Validations in-run: M11 against run 11's
     kernel-domain block; g, u_1..u_12 against Heart's arb values;
     diag excess D_j = M_jj - g w_j^2 and lambda_min(M22), Schur
     sweep against Heart's float sweep (GPPDiscovery5 2d3506b).
  3. The tail loss, in BOTH pictures (plain Schur S(K); bordered
     Schur SB(K) of S~ = M - uu^T/g): Lambda(60) measured as
     lambda_min SB(60) - lambda_min SB(800) (identical quadrature
     throughout, so systematics cancel in the difference); the
     beyond-812 tail bounded from the measured increment ratios;
     the residual coupling v~ after border+Y_72 elimination measured
     column by column: C_v = max_j ||v~_.j||_2 w_j^3, its rank-one
     p-profile part q (p_j = w_j/w_j^2) against the prediction
     q = -d + (h(l)/g) u, and the leftover constant C_res.
  4. A second quadrature resolution (64 nodes/period) recomputes
     lambda_min S(60), S(800), SB(60), SB(800): the drift of
     Lambda(60) is the empirical error bar.

No zero data anywhere.  Epistemic status: NUMERICAL (float64,
independent code path), with the analytic steps labelled in the
notes.  This script certifies a DIFFERENCE of Schur complements at
L = 0.5 only; it is not a proof of OddFace(0.5) (that needs (a)
plus this bound plus the assembly), and it says nothing about
other L.
"""

import json
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import digamma

T_START = time.time()
OUT = {'label': 'NUMERICAL (float64, independent code path)'}
ELL = 0.25
SQ2L = np.sqrt(2.0 / ELL)                      # sqrt(8)
BIGPSI_I2 = 16.0 * np.sinh(0.25) ** 2          # Psi(i/2) = 1.02101...
LAM = {2: np.log(2), 3: np.log(3), 4: np.log(2), 5: np.log(5),
       7: np.log(7), 8: np.log(2), 9: np.log(3), 11: np.log(11)}
JMAX = 812
K_OF = {20: 20, 50: 50, 60: 60, 100: 100, 200: 200, 400: 400, 800: 800}

# Heart's filed values (GPPDiscovery5 2d3506b, p2_schur_sweep_results
# .json and p2_bordered_arb_results.json; 3ee15a8 for P1) -- comparison
# targets only, never inputs.
HEART = {
    'g': 1.23787e-4,
    'u12': [-3.443378638e-5, -3.052769456e-5, 9.68767956e-6,
            -8.10263991e-6, 6.56443284e-6, -5.28510625e-6,
            4.49691654e-6, -3.7203197e-6, 3.2564694e-6, -2.9327152e-6,
            2.6454502e-6, -2.4193589e-6],
    'beta0': 1.30163e-8,
    'leading_limit': 1.4462598791127561e-08,
    'sweep': {  # K: (lam_min_M22, lam_min_schur, diag_excess_first, diag_excess_last)
        20: (4.671221593413881e-10, 1.4710204512301735e-08,
             1.4890079691675708e-08, 5.184939761522865e-10),
        50: (4.098055573924145e-11, 1.4340862021910105e-08,
             1.4890079691675708e-08, 4.9089877095283486e-11),
        100: (3.951323671534351e-12, 1.4233979154158182e-08,
              1.4890079691675708e-08, 4.806264001073961e-12),
        200: (3.202053409787337e-13, 1.4189879276372636e-08,
              1.4890079691675708e-08, 4.227772184840721e-13),
        400: (2.3535143585957793e-14, 1.417125506096058e-08,
              1.4890079691675708e-08, 3.2208545267411224e-14)},
    'archpp0': 1.52192,          # archTerm''(0), OF1 (Round 4)
    'arch0': -1.54293308589872,  # banked accurate archTerm(0)
    'pole0': 2.0849135341,       # poleTerm(0) = 2*bigPsiI2^2
}

# ============================================================ Part C-lite
# arch values/derivatives by period-wise Gauss-Legendre (run-11 rig).
def arch_integral(T, nodes_per, dgrid, power, kind):
    """(1/pi) int_0^T Psi^2 w(t) t^power * trig(t d) dt on a d-grid.
    kind = 'cos' (even derivatives) or 'sin' (odd)."""
    gx, gw = leggauss(nodes_per)
    nper = int(np.ceil(T / (2 * np.pi)))
    edges = 2 * np.pi * np.arange(nper + 1)
    cen = 0.5 * (edges[:-1] + edges[1:])
    half = 0.5 * (edges[1:] - edges[:-1])
    tn = (cen[:, None] + half[:, None] * gx[None, :]).ravel()
    tw = (half[:, None] * gw[None, :]).ravel()
    psi2 = (np.sin(tn / 2) / (tn / 2)) ** 4
    wgt = digamma(0.25 + 0.5j * tn).real - np.log(np.pi)
    base = tw * psi2 * wgt * tn ** power / np.pi
    out = np.empty(len(dgrid))
    for i0 in range(0, len(dgrid), 200):
        dg = dgrid[i0:i0 + 200]
        trig = np.cos(np.outer(tn, dg)) if kind == 'cos' \
            else np.sin(np.outer(tn, dg))
        out[i0:i0 + 200] = base @ trig
    return out


def triG(y):
    a = np.abs(y)
    return np.where(a <= 1, 2 / 3 - y ** 2 + a ** 3 / 2,
                    np.where(a <= 2, (2 - a) ** 3 / 6, 0.0))


def triG2(y):
    """triG''(y): -2 + 3|y| (|y|<=1); (2-|y|) (1<=|y|<=2); 0 else."""
    a = np.abs(y)
    return np.where(a <= 1, -2 + 3 * a,
                    np.where(a <= 2, 2 - a, 0.0))


def triG3(y):
    """triG'''(y): 3 sign(y) (|y|<1); -sign(y) (1<|y|<2); 0 else."""
    a = np.abs(y)
    s = np.sign(y)
    return np.where(a < 1, 3 * s, np.where(a < 2, -s, 0.0))


def prime_S(d, deriv=0):
    d = np.asarray(d, dtype=float)
    tot = np.zeros_like(d)
    for n, lam in LAM.items():
        ln = np.log(n)
        if deriv == 0:
            f = triG(ln + d) + triG(-ln + d)
        elif deriv == 2:
            f = triG2(ln + d) + triG2(-ln + d)
        elif deriv == 3:
            f = triG3(ln + d) + triG3(-ln + d)
        tot = tot + lam / np.sqrt(n) * f
    return tot


B2 = BIGPSI_I2 ** 2
OUT['pole0_check'] = {'mine': float(2 * B2), 'banked': HEART['pole0'],
                      'absdiff': float(abs(2 * B2 - HEART['pole0']))}

# arch(0), arch(0.5) at the run-11 resolution for g.
_d2 = np.array([0.0, 0.5])
arch_d2 = arch_integral(4000.0, 64, _d2, 0, 'cos')
OUT['arch0_mine'] = float(arch_d2[0])
OUT['arch0_absdiff'] = float(abs(arch_d2[0] - HEART['arch0']))
KS0 = 2 * B2 - float(prime_S(np.array([0.0]))[0]) + arch_d2[0]
KS05 = 2 * B2 * np.cosh(0.25) - float(prime_S(np.array([0.5]))[0]) \
    + arch_d2[1]
G_VAL = 0.5 * (KS0 - KS05)
OUT['g_kernel'] = float(G_VAL)
OUT['g_heart'] = HEART['g']
OUT['g_reldiff'] = float(abs(G_VAL - HEART['g']) / HEART['g'])
print('g (kernel):', G_VAL, ' Heart:', HEART['g'],
      ' rel diff', OUT['g_reldiff'])
print('archTerm(0): mine', arch_d2[0], ' banked', HEART['arch0'],
      ' |diff|', OUT['arch0_absdiff'])

# K_S''(d) on [0, 0.5]: pole'' exact, prime'' exact, arch'' by GL.
DGRID2 = np.linspace(0.0, 0.5, 801)
archpp = arch_integral(20000.0, 40, DGRID2, 2, 'cos')
archpp = -archpp                                   # d^2 cos -> -t^2
polepp = 0.5 * B2 * np.cosh(DGRID2 / 2)
KPP = polepp - prime_S(DGRID2, deriv=2) + archpp
OUT['archpp0_mine'] = float(archpp[0])
OUT['archpp0_heart'] = HEART['archpp0']
OUT['archpp0_absdiff'] = float(abs(archpp[0] - HEART['archpp0']))
print("arch''(0): mine", archpp[0], ' Heart OF1:', HEART['archpp0'],
      ' |diff|', OUT['archpp0_absdiff'],
      ' (tail allowance 3.8e-3)')
KPP0 = float(KPP[0])
KPP05 = float(KPP[-1])
H_ELL = 0.5 * (KPP0 - KPP05)                       # b2(l)
OUT['Kpp0'] = KPP0
OUT['Kpp05'] = KPP05
OUT['h_ell'] = float(H_ELL)
print("K''(0) =", KPP0, " K''(0.5) =", KPP05, ' h(l) = b2(l) =', H_ELL)

# b2(u) = (1/2)[K''(l-u) - K''(l+u)] on a fine u-grid; sine transforms.
NU = 32768
ugrid = np.linspace(0.0, ELL, NU + 1)
KPP_interp = lambda x: np.interp(x, DGRID2, KPP)
b2_u = 0.5 * (KPP_interp(ELL - ugrid) - KPP_interp(ELL + ugrid))
du = ugrid[1] - ugrid[0]
trap = np.full(NU + 1, du); trap[0] = du / 2; trap[-1] = du / 2
OUT['b2_sup'] = float(np.max(np.abs(b2_u)))
OUT['b2_L1'] = float(np.sum(trap * np.abs(b2_u)))
OUT['b2_L2'] = float(np.sqrt(np.sum(trap * b2_u ** 2)))
print('|b2| sup/L1/L2:', OUT['b2_sup'], OUT['b2_L1'], OUT['b2_L2'])

om_all = 4 * np.pi * np.arange(1, JMAX + 1)
w_all = SQ2L * (-1.0) ** (np.arange(1, JMAX + 1) + 1) / om_all
# sine transforms of b2: shat_j = int_0^l sin(w_j u) b2(u) du
shat = np.empty(JMAX)
for i0 in range(0, JMAX, 128):
    jj = np.arange(i0 + 1, min(i0 + 128, JMAX) + 1)
    S = np.sin(np.outer(om_all[i0:i0 + len(jj)], ugrid))
    shat[i0:i0 + len(jj)] = (S * trap[None, :]) @ b2_u
u_all = G_VAL * w_all - SQ2L * shat / om_all ** 2
OUT['u_spot'] = {f'u{j+1}': {'mine': float(u_all[j]),
                             'heart_arb': HEART['u12'][j],
                             'absdiff': float(abs(u_all[j]
                                                  - HEART['u12'][j]))}
                 for j in range(12)}
print('u_1..u_12 max |diff| vs Heart arb:',
      max(v['absdiff'] for v in OUT['u_spot'].values()))
# d_i = <phi_i, b2>, i <= 12 (same transforms, phi_i = sqrt(2/l) sin)
d_vec = SQ2L * shat[:12]
OUT['d_vec'] = [float(x) for x in d_vec]
print('d_i (i<=12):', ['%.4e' % x for x in d_vec])

# C3 = int_0^l || d^3_s G(s,.) ||_2 ds from K'''(d), d in [-0.5, 0.5].
DGRID3 = np.linspace(0.0, 0.5, 2001)
archppp = arch_integral(20000.0, 32, DGRID3, 3, 'sin')
poleppp = 0.25 * B2 * np.sinh(DGRID3 / 2)
KPPP_half = poleppp - prime_S(DGRID3, deriv=3) + archppp   # d >= 0
# log-model fit for arch''' on d in [0.004, 0.04] (arch''' ~ A log d + B
# there; the fit feeds the singular-cell integrals of C3).
fit_m = (DGRID3 >= 0.004) & (DGRID3 <= 0.04)
Afit, Bfit = np.polyfit(np.log(DGRID3[fit_m]), archppp[fit_m], 1)
OUT['archppp_logfit'] = {'A': float(Afit), 'B': float(Bfit),
                         'max_resid_fit_range':
                         float(np.max(np.abs(
                             archppp[fit_m]
                             - (Afit * np.log(DGRID3[fit_m]) + Bfit))))}
print('archppp log fit: A =', Afit, ' B =', Bfit,
      ' max resid', OUT['archppp_logfit']['max_resid_fit_range'])


def KPPP(d):
    """K_S'''(d) for arbitrary d (odd extension; log model |d|<0.004
    applied to the arch part only via the full-value fit)."""
    d = np.asarray(d, dtype=float)
    sgn = np.sign(d)
    a = np.abs(d)
    val = np.interp(a, DGRID3, KPPP_half)
    small = a < 0.004
    if np.any(small):
        # replace arch part by fit: archppp_fit - interp'd archppp
        arch_i = np.interp(a[small], DGRID3, archppp)
        arch_f = Afit * np.log(np.maximum(a[small], 1e-12)) + Bfit
        val[small] = val[small] - arch_i + arch_f
    return sgn * val


sgrid = np.linspace(0.0, ELL, 51)
uu = np.linspace(0.0, ELL, 801)
duu = uu[1] - uu[0]
tw_u = np.full(len(uu), duu); tw_u[0] = duu / 2; tw_u[-1] = duu / 2
norms = np.empty(len(sgrid))
for i, s in enumerate(sgrid):
    diff = KPPP(s - uu) - KPPP(s + uu)
    norms[i] = 0.5 * np.sqrt(np.sum(tw_u * diff ** 2))
C3 = float(np.trapezoid(norms, sgrid))
OUT['C3'] = C3
C_R = SQ2L * (OUT['b2_L2'] + C3)
OUT['C_R'] = float(C_R)
print('C3 =', C3, ' C_R = sqrt(2/l)(||b2||_2 + C3) =', C_R)
print('Part C elapsed %.1fs' % (time.time() - T_START))


# ============================================================ Part A
# Full Gram M on modes 1..812 by period-wise Gauss-Legendre.
# (First assembly run had S_j too small by exactly l = 1/4 -- an extra
# ELL factor against sinc(x) = sin x / x; caught by the S_j direct-
# integral check in OUT['S_checks'], fixed, rerun.  Kernel-side Part C
# was unaffected.  Final run: S checks ~1e-16, M11 vs run-11 kernel
# matrix max |entry diff| 1.6e-12.)
def assemble_M(nodes_per, T_extra=400.0):
    gx, gw = leggauss(nodes_per)
    T = om_all[-1] + T_extra
    nper = int(np.ceil(T / (2 * np.pi)))
    edges = 2 * np.pi * np.arange(nper + 1)
    cen = 0.5 * (edges[:-1] + edges[1:])
    half = 0.5 * (edges[1:] - edges[:-1])
    tn = (cen[:, None] + half[:, None] * gx[None, :]).ravel()
    tw = (half[:, None] * gw[None, :]).ravel()
    psi2 = (np.sin(tn / 2) / (tn / 2)) ** 4
    wgt = digamma(0.25 + 0.5j * tn).real - np.log(np.pi)
    Pt = np.zeros_like(tn)
    for n, lam in LAM.items():
        Pt += lam / np.sqrt(n) * np.cos(tn * np.log(n))
    V = wgt - 2 * Pt
    Wq = tw * psi2 * V / np.pi
    sig = -(-1.0) ** np.arange(1, JMAX + 1) * om_all \
        * np.sinh(ELL / 2) / (0.25 + om_all ** 2)
    M = np.empty((JMAX, JMAX))
    CH = 80
    for i0 in range(0, JMAX, CH):
        i1 = min(i0 + CH, JMAX)
        omr = om_all[i0:i1]
        Sr = omr[:, None] * (ELL * np.sinc(
            (tn[None, :] - omr[:, None]) * ELL / np.pi)) \
            / (tn[None, :] + omr[:, None])
        SrW = Sr * Wq[None, :]
        for j0 in range(i0, JMAX, CH):
            j1 = min(j0 + CH, JMAX)
            omc = om_all[j0:j1]
            Sc = omc[:, None] * (ELL * np.sinc(
                (tn[None, :] - omc[:, None]) * ELL / np.pi)) \
                / (tn[None, :] + omc[:, None])
            blk = (2 / ELL) * (SrW @ Sc.T)
            blk -= (2 / ELL) * 2 * B2 * np.outer(sig[i0:i1], sig[j0:j1])
            M[i0:i1, j0:j1] = blk
            if j0 > i0:
                M[j0:j1, i0:i1] = blk.T
    return 0.5 * (M + M.T), sig


print('assembling M (48 nodes/period) ...')
M, sig_all = assemble_M(48)
print('M assembled, elapsed %.1fs' % (time.time() - T_START))
np.save('/tmp/flint12_M.npy', M)

# closed-form spot validations: sigma_j and S_j against direct GL.
gxq, gwq = leggauss(200)
xn = 0.5 * ELL * (gxq + 1)
xw = 0.5 * ELL * gwq
sig_check = {}
for j in (1, 2, 13):
    om = 4 * np.pi * j
    num = np.sum(xw * np.sinh(xn / 2) * np.sin(om * xn))
    sig_check[f'sigma{j}'] = {'direct': float(num),
                              'closed': float(sig_all[j - 1]),
                              'absdiff': float(abs(num - sig_all[j - 1]))}
OUT['sigma_checks'] = sig_check
S_check = {}
for (j, t) in ((3, 40.0), (13, 170.0), (100, 1300.0)):
    om = 4 * np.pi * j
    num = np.sum(xw * np.sin(om * xn) * np.sin(t * xn))
    closed = om * (ELL * np.sinc((t - om) * ELL / np.pi)) / (t + om)
    S_check[f'S{j}_at_{t}'] = {'direct': float(num),
                               'closed': float(closed),
                               'absdiff': float(abs(num - closed))}
OUT['S_checks'] = S_check
print('sigma/S closed-form checks done')

# ============================================================ Part B
M11 = M[:12, :12]
ev11 = np.linalg.eigvalsh(M11)
OUT['lambda_min_M11_mine'] = float(ev11[0])
r11 = json.load(open('/home/hatch/workspace/gppdiscovery3/research/'
                     'flint_p1_kernel_check_2026_11_15_results.json'))
M11_prev = np.array(r11['M_mine'])
OUT['M11_vs_run11'] = {
    'max_absdiff': float(np.max(np.abs(M11 - M11_prev))),
    'lambda_min_run11': r11['lambda_min_mine'],
    'lambda_min_absdiff': float(abs(ev11[0] - r11['lambda_min_mine']))}
print('M11 lambda_min mine:', ev11[0], ' run-11 kernel:',
      r11['lambda_min_mine'], ' max entry |diff| vs run 11:',
      OUT['M11_vs_run11']['max_absdiff'])

wvec = w_all
D_ex = np.array([M[j, j] - G_VAL * wvec[j] ** 2
                 for j in range(12, JMAX)])
OUT['diag_excess'] = {
    'D13': float(D_ex[0]),
    'D13_heart': HEART['sweep'][400][2],
    'D412': float(D_ex[399]),
    'D412_heart': HEART['sweep'][400][3],
    'D812': float(D_ex[-1]),
    'D_times_om4_over_logom_first': float(
        D_ex[0] * om_all[12] ** 4 / np.log(om_all[12])),
    'D_times_om4_over_logom_last': float(
        D_ex[-1] * om_all[-1] ** 4 / np.log(om_all[-1]))}
print('diag excess: D13 mine', D_ex[0], ' Heart', HEART['sweep'][400][2],
      '; D412 mine', D_ex[399], ' Heart', HEART['sweep'][400][3])


def schur_lambda_min(Mfull, K):
    A = Mfull[:12, 12:12 + K]
    KK = Mfull[12:12 + K, 12:12 + K]
    S = Mfull[:12, :12] - A @ np.linalg.solve(KK, A.T)
    return float(np.linalg.eigvalsh(0.5 * (S + S.T))[0]), S


sweep = {}
for K in (20, 50, 60, 100, 200, 400, 800):
    lam22 = float(np.linalg.eigvalsh(M[12:12 + K, 12:12 + K])[0])
    lamS, _ = schur_lambda_min(M, K)
    row = {'lambda_min_M22': lam22, 'lambda_min_schur': lamS}
    if K in HEART['sweep']:
        row['heart_lambda_min_M22'] = HEART['sweep'][K][0]
        row['heart_lambda_min_schur'] = HEART['sweep'][K][1]
        row['schur_absdiff_vs_heart'] = abs(lamS - HEART['sweep'][K][1])
    sweep[str(K)] = row
    print(f'K={K}: lam M22 {lam22:.6e}  lam Schur {lamS:.10e}',
          (' Heart schur %.10e' % HEART['sweep'][K][1]) if K in
          HEART['sweep'] else '')
OUT['sweep_plain'] = sweep

# ============================================================ Part D
# Bordered picture: S~ = M - u u^T / g on all modes; stage-K bordered
# Schur = plain stage-K Schur of S~.
u_vec = u_all
St = M - np.outer(u_vec, u_vec) / G_VAL
lead = M11 - np.outer(u_vec[:12], u_vec[:12]) / G_VAL
OUT['leading_limit_mine'] = float(np.linalg.eigvalsh(lead)[0])
OUT['leading_limit_heart'] = HEART['leading_limit']
print('bordered leading limit: mine', OUT['leading_limit_mine'],
      ' Heart float', HEART['leading_limit'])

sweep_b = {}
for K in (20, 50, 60, 100, 200, 400, 800):
    lamSB, _ = schur_lambda_min(St, K)
    sweep_b[str(K)] = {'lambda_min_schur_bordered': lamSB}
    print(f'K={K}: lam Schur bordered {lamSB:.10e}')
OUT['sweep_bordered'] = sweep_b

LAM60 = sweep_b['60']['lambda_min_schur_bordered'] \
    - sweep_b['800']['lambda_min_schur_bordered']
LAM60_plain = sweep['60']['lambda_min_schur'] \
    - sweep['800']['lambda_min_schur']
OUT['Lambda60_bordered_measured'] = float(LAM60)
OUT['Lambda60_plain_measured'] = float(LAM60_plain)
print('Lambda(60) measured: bordered', LAM60, ' plain', LAM60_plain)

# increment ratios of the bordered sweep (tail extrapolation)
incr = {}
Ks = [50, 100, 200, 400]
for a, b in zip(Ks[:-1], Ks[1:]):
    incr[f'{a}->{b}'] = sweep_b[str(a)]['lambda_min_schur_bordered'] \
        - sweep_b[str(b)]['lambda_min_schur_bordered']
ratios = [incr['100->200'] / incr['50->100'],
          incr['200->400'] / incr['100->200']]
r_use = max(ratios) * 1.25
tail800 = incr['200->400'] * r_use / (1 - r_use)   # bound for >800 loss
OUT['increments_bordered'] = {k: float(v) for k, v in incr.items()}
OUT['increment_ratios'] = [float(x) for x in ratios]
OUT['tail_beyond_800_bound'] = float(tail800)
OUT['Lambda60_total_bound'] = float(LAM60 + tail800)
OUT['beta0'] = HEART['beta0']
OUT['room_factor_vs_beta0'] = float(HEART['beta0']
                                    / (LAM60 + tail800))
print('increments:', incr, ' ratios:', ratios)
print('tail beyond 800 (ratio bound):', tail800)
print('Lambda(60) total bound:', LAM60 + tail800,
      ' vs beta0', HEART['beta0'], ' room x',
      OUT['room_factor_vs_beta0'])

# residual coupling v~ after border + Y_72 elimination; Z = modes > 72.
Y = slice(12, 72)
Z = slice(72, JMAX)
S11 = St[:12, :12]
SY = St[Y, Y]
v_tilde = St[:12, Z] - St[:12, Y] @ np.linalg.solve(SY, St[Y, Z])
SZ_eff = St[Z, Z] - St[Z, Y] @ np.linalg.solve(SY, St[Y, Z])
SZ_eff = 0.5 * (SZ_eff + SZ_eff.T)
colnorm = np.linalg.norm(v_tilde, axis=0)
omZ = om_all[72:]
Cv = float(np.max(colnorm * omZ ** 3))
jworst = int(np.argmax(colnorm * omZ ** 3)) + 73
OUT['C_vtilde'] = {'C': Cv, 'worst_mode': jworst,
                   'colnorm_at_73': float(colnorm[0]),
                   'colnorm_at_412': float(colnorm[399]),
                   'colnorm_at_812': float(colnorm[-1])}
# rank-one p-profile part: p_j = w_j / w_j^2 over Z
pZ = wvec[72:] / omZ ** 2
q_hat = (v_tilde @ pZ) / (pZ @ pZ)
resid = v_tilde - np.outer(q_hat, pZ)
C_res = float(np.max(np.linalg.norm(resid, axis=0) * omZ ** 3))
q_pred = -d_vec + (H_ELL / G_VAL) * u_vec[:12]
OUT['p_profile'] = {
    'q_hat': [float(x) for x in q_hat],
    'q_pred_minus_d_plus_h_over_g_u': [float(x) for x in q_pred],
    'q_hat_norm': float(np.linalg.norm(q_hat)),
    'q_pred_norm': float(np.linalg.norm(q_pred)),
    'cosine_qhat_qpred': float(q_hat @ q_pred
                               / (np.linalg.norm(q_hat)
                                  * np.linalg.norm(q_pred))),
    'C_res_after_p_profile': C_res}
print('v~: C_v =', Cv, 'at mode', jworst,
      '; after p-profile removal C_res =', C_res)
print('q_hat norm', np.linalg.norm(q_hat), ' q_pred norm',
      np.linalg.norm(q_pred), ' cosine',
      OUT['p_profile']['cosine_qhat_qpred'])
# loss quadratic form and floors
loss_mat = v_tilde @ np.linalg.solve(SZ_eff, v_tilde.T)
loss_mat = 0.5 * (loss_mat + loss_mat.T)
OUT['loss_72_to_812'] = {
    'op_norm': float(np.linalg.eigvalsh(loss_mat)[-1]),
    'lambda_min_effect_check': float(
        sweep_b['60']['lambda_min_schur_bordered']
        - np.linalg.eigvalsh(
            (S11 - St[:12, Y] @ np.linalg.solve(SY, St[Y, :12]))
            - loss_mat)[0]),
    'direct_lambda_min_full_bordered': sweep_b['800']
        ['lambda_min_schur_bordered']}
OUT['floors'] = {
    'lambda_min_SZeff_73_812': float(np.linalg.eigvalsh(SZ_eff)[0]),
    'lambda_min_SZeff_413_812_block': float(np.linalg.eigvalsh(
        SZ_eff[340:, 340:])[0]),
    'lambda_min_St_413_812': float(np.linalg.eigvalsh(St[412:, 412:])[0])}
print('loss 72->812 op norm:', OUT['loss_72_to_812']['op_norm'],
      ' floors:', OUT['floors'])

# analytic tail beyond 812 from the measured constants:
# ||v~_.j|| <= C_v w_j^-3 exactly as measured shape; sum w_j^-6.
sum_om6_tail = float(np.sum(1.0 / om_all ** 6)
                     + 1.0 / (5 * om_all[-1] ** 5 * 4 * np.pi))
# (sum over j<=812 computed; tail by integral bound (4pi)^-6/(5*812^5)
#  = 1/(5 * om_812^5 * 4pi) * (4pi)^-1... kept explicit:)
sum_om6_tail = float((4 * np.pi) ** -6
                     * (1.0 / (5 * 812.0 ** 5) + 1.0 / 812.0 ** 6))
OUT['tail_beyond_812_analytic'] = {
    'sum_omega^-6_beyond_812': sum_om6_tail,
    'bound_Cv2_times_sum_over_floor': float(
        12 * Cv ** 2 * sum_om6_tail
        / OUT['floors']['lambda_min_SZeff_413_812_block'])}
print('analytic >812 tail (Cv^2 * sum / floor):',
      OUT['tail_beyond_812_analytic']['bound_Cv2_times_sum_over_floor'])

# ============================================================ Part E
# Second quadrature resolution (64 nodes/period): drift of the
# headline differences.  u, g are kernel-side (resolution-free here).
print('assembling M at 64 nodes/period (drift check) ...')
M2, _ = assemble_M(64)
St2 = M2 - np.outer(u_vec, u_vec) / G_VAL
drift = {}
for K in (60, 800):
    lS2, _ = schur_lambda_min(M2, K)
    lSB2, _ = schur_lambda_min(St2, K)
    drift[f'K{K}'] = {'plain_res64': lS2, 'bordered_res64': lSB2}
LAM60_b_res64 = drift['K60']['bordered_res64'] \
    - drift['K800']['bordered_res64']
LAM60_p_res64 = drift['K60']['plain_res64'] - drift['K800']['plain_res64']
OUT['res64_drift'] = {
    **drift,
    'Lambda60_bordered_res64': float(LAM60_b_res64),
    'Lambda60_plain_res64': float(LAM60_p_res64),
    'Lambda60_bordered_drift': float(abs(LAM60_b_res64 - LAM60)),
    'Lambda60_plain_drift': float(abs(LAM60_p_res64 - LAM60_plain))}
print('res64: Lambda60 bordered', LAM60_b_res64, ' drift',
      OUT['res64_drift']['Lambda60_bordered_drift'],
      '; plain', LAM60_p_res64, ' drift',
      OUT['res64_drift']['Lambda60_plain_drift'])

OUT['runtime_seconds'] = time.time() - T_START
with open('/home/hatch/workspace/gppdiscovery3/research/'
          'flint_tail_bound_2026_11_15_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('done', OUT['runtime_seconds'])
