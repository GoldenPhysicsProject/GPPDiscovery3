#!/usr/bin/env python
"""Flint thirteenth run (2026-11-15): the reached-subspace floor for
S_Z^eff -- the proof-grade step named in run 12 section 2 Step 4.

CONTEXT.  Run 12 (flint_tail_bound_2026_11_15) certified the OddFace
(0.5) tail bound (b) NUMERICALLY: Lambda(60) <= 4.43e-11 against
beta0 = 1.30163e-8.  Its closed analytic route FAILED, located
exactly: the route charges the tail loss through the global floor
lambda_min(S_Z^eff) = 1.79e-15, which is D_812-scale (attained at
the near-null TOP modes of the truncation), while the coupling
reaches those modes only at noise-floor strength.  The measured
loss operator L = v~ (S_Z^eff)^{-1} v~^T has lambda_max = 5.87e-8
but its action on the bottom eigendirection of the stage-60
bordered form -- the only action the assembly consumes -- is
1.62e-11.  Run 12 Step 4 named the missing proof step: a floor
for S_Z^eff restricted to the coupling-reached subspace
range(v~^T), worth ~3.6e3 numerically.  THIS script computes that
floor and the thresholds it must beat.

OBJECTS (run 12's assembly, self-contained).  L = 0.5, l = 1/4,
sine basis phi_j, omega_j = 4 pi j.  M = full Gram (spectral,
period-wise Gauss-Legendre, 48 nodes/period, modes 1..812);
S~ = M - u u^T/g the bordered form (g = G(l,l), u_j = B(phi_j,
delta_l), kernel side).  Partition X = modes 1..12, Y = modes
13..72, Z = modes 73..812.  After border + Y-elimination:
    v~     = S~_XZ - S~_XY S~_YY^{-1} S~_YZ   (12 x |Z| coupling)
    A      = S_Z^eff = S~_ZZ - S~_ZY S~_YY^{-1} S~_YZ  (Z-block)
    A60    = S~_XX - S~_XY S~_YY^{-1} S~_YX   (stage-60 form)
    L      = v~ A^{-1} v~^T                    (loss form, 12 x 12)
    SB(60) = lambda_min(A60),  v0 its bottom eigenvector,
    Lambda(60) = SB(60) - lambda_min(A60 - L)  (exact identity).

THRESHOLDS (derived before the floors, in the output order):
    SigmaF  = ||v~||_F^2 = sum_{j in Z} ||v~_.j||^2,
    SigmaV0 = ||v~^T v0||^2,
    c*_F    = SigmaF / beta0   (floor making SigmaF/c < beta0),
    c*_v0   = SigmaV0 / beta0  (floor making SigmaV0/c < beta0),
    c*_tail = 12 C_v^2 sum_{j>812} omega_j^-6 / beta0
              (run 12's filed crude-route threshold, = 6.12e-15).

FLOORS (at truncations J = 412, 612, 812; Z_J = modes 73..J):
    lam_min(A_J)                    global floor,
    lam_min(U_J^T A_J U_J)          compression floor,
    c_inv(J) = lam_min((U_J^T A_J^{-1} U_J)^{-1})
                                    inverse-form (Schur) floor,
    c_dir(J) = SigmaV0 / (v0^T L_J v0)   directional floor,
with U_J an orthonormal basis of range(v~_J^T) (SVD).  The loss
route: v0^T L v0 = <v~^T v0, A^{-1} v~^T v0> <= SigmaV0 / c_inv
<= SigmaF / c_inv (Rayleigh bound on the reached subspace).

Also: near-null eigendirections of A (bottom 8) with band-mass
profiles and overlaps with range(v~^T); the reached profile of
y = v~^T v0; D-normalised off-diagonal size rho_Z of A; the
diagonal profile A_jj * omega_j^4 / log omega_j; the Y-elimination
correction to the Z diagonal; SigmaD = sum_j ||v~_.j||^2 / A_jj;
lambda_min(v~ v~^T); the analytic chain bound
c_inv >= (1 - rho_Z) * lambda_min(v~ v~^T) / SigmaD; the
generalised eigenvalue mu_gen = lam_max(A60^{-1/2} L A60^{-1/2})
with Lambda(60) <= mu_gen * SB(60); and a 64-nodes/period drift
check of the J = 812 floors.

No zero data anywhere.  Epistemic status: NUMERICAL (float64),
same code path as run 12, validated in-run against run 12's filed
values (lambda_min M11, SB(60), SB(800), lambda_max(L),
lambda_min(S_Z^eff), Lambda(60) via L).
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
BETA0 = 1.30163e-8        # Heart R6 leading-order certified floor
BETA1 = 1.41006e-8        # Heart (a) at K = 60, audited (Turbo relay)

HEART = {
    'g': 1.23787e-4,
    'u12': [-3.443378638e-5, -3.052769456e-5, 9.68767956e-6,
            -8.10263991e-6, 6.56443284e-6, -5.28510625e-6,
            4.49691654e-6, -3.7203197e-6, 3.2564694e-6, -2.9327152e-6,
            2.6454502e-6, -2.4193589e-6],
    'arch0': -1.54293308589872,
    'pole0': 2.0849135341}

# ============================================================ Part C-lite
# (verbatim from run 12: kernel-side g, u_all, d_vec)
def arch_integral(T, nodes_per, dgrid, power, kind):
    """(1/pi) int_0^T Psi^2 w(t) t^power * trig(t d) dt on a d-grid."""
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
    a = np.abs(y)
    return np.where(a <= 1, -2 + 3 * a,
                    np.where(a <= 2, 2 - a, 0.0))


def prime_S(d, deriv=0):
    d = np.asarray(d, dtype=float)
    tot = np.zeros_like(d)
    for n, lam in LAM.items():
        ln = np.log(n)
        if deriv == 0:
            f = triG(ln + d) + triG(-ln + d)
        elif deriv == 2:
            f = triG2(ln + d) + triG2(-ln + d)
        tot = tot + lam / np.sqrt(n) * f
    return tot


B2 = BIGPSI_I2 ** 2
_d2 = np.array([0.0, 0.5])
arch_d2 = arch_integral(4000.0, 64, _d2, 0, 'cos')
KS0 = 2 * B2 - float(prime_S(np.array([0.0]))[0]) + arch_d2[0]
KS05 = 2 * B2 * np.cosh(0.25) - float(prime_S(np.array([0.5]))[0]) \
    + arch_d2[1]
G_VAL = 0.5 * (KS0 - KS05)
OUT['g_kernel'] = float(G_VAL)
OUT['g_reldiff_vs_heart'] = float(abs(G_VAL - HEART['g']) / HEART['g'])
print('g (kernel):', G_VAL, ' rel diff vs Heart', OUT['g_reldiff_vs_heart'])

DGRID2 = np.linspace(0.0, 0.5, 801)
archpp = arch_integral(20000.0, 40, DGRID2, 2, 'cos')
archpp = -archpp
polepp = 0.5 * B2 * np.cosh(DGRID2 / 2)
KPP = polepp - prime_S(DGRID2, deriv=2) + archpp
KPP0 = float(KPP[0])
KPP05 = float(KPP[-1])
H_ELL = 0.5 * (KPP0 - KPP05)
OUT['h_ell'] = float(H_ELL)

NU = 32768
ugrid = np.linspace(0.0, ELL, NU + 1)
KPP_interp = lambda x: np.interp(x, DGRID2, KPP)
b2_u = 0.5 * (KPP_interp(ELL - ugrid) - KPP_interp(ELL + ugrid))
du = ugrid[1] - ugrid[0]
trap = np.full(NU + 1, du); trap[0] = du / 2; trap[-1] = du / 2

om_all = 4 * np.pi * np.arange(1, JMAX + 1)
w_all = SQ2L * (-1.0) ** (np.arange(1, JMAX + 1) + 1) / om_all
shat = np.empty(JMAX)
for i0 in range(0, JMAX, 128):
    jj = np.arange(i0 + 1, min(i0 + 128, JMAX) + 1)
    S = np.sin(np.outer(om_all[i0:i0 + len(jj)], ugrid))
    shat[i0:i0 + len(jj)] = (S * trap[None, :]) @ b2_u
u_all = G_VAL * w_all - SQ2L * shat / om_all ** 2
d_vec = SQ2L * shat[:12]
OUT['u12_max_absdiff_vs_heart'] = float(
    np.max(np.abs(u_all[:12] - np.array(HEART['u12']))))
print('u_1..u_12 max |diff| vs Heart arb:', OUT['u12_max_absdiff_vs_heart'])
print('Part C elapsed %.1fs' % (time.time() - T_START))

# ============================================================ Part A
# Full Gram M on modes 1..812 (verbatim from run 12).
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

# ============================================================ Part F
# The reached-subspace floor.  X = modes 1..12, Y = 13..72, Z = 73..812
# (0-based slices as in run 12).
R12 = json.load(open('/home/hatch/workspace/gppdiscovery3/research/'
                     'flint_tail_bound_2026_11_15_results.json'))
St = M - np.outer(u_all, u_all) / G_VAL
X = slice(0, 12)
Y = slice(12, 72)
Z = slice(72, JMAX)
SY = St[Y, Y]
SYinv_StYZ = np.linalg.solve(SY, St[Y, Z])
SYinv_StYX = np.linalg.solve(SY, St[Y, X])
v_tilde = St[X, Z] - St[X, Y] @ SYinv_StYZ          # 12 x 740
SZ_eff = St[Z, Z] - St[Z, Y] @ SYinv_StYZ
SZ_eff = 0.5 * (SZ_eff + SZ_eff.T)
A60 = St[X, X] - St[X, Y] @ SYinv_StYX
A60 = 0.5 * (A60 + A60.T)

ev60, evec60 = np.linalg.eigh(A60)
SB60 = float(ev60[0])
v0 = evec60[:, 0]
gap60 = float(ev60[1] - ev60[0])
loss_mat = v_tilde @ np.linalg.solve(SZ_eff, v_tilde.T)
loss_mat = 0.5 * (loss_mat + loss_mat.T)
lam_max_L = float(np.linalg.eigvalsh(loss_mat)[-1])
v0Lv0 = float(v0 @ loss_mat @ v0)
full_effect = SB60 - float(np.linalg.eigvalsh(A60 - loss_mat)[0])
lam_min_A = float(np.linalg.eigvalsh(SZ_eff)[0])
resid_Lv0 = float(np.linalg.norm(loss_mat @ v0 - v0Lv0 * v0))
# generalised eigenvalue mu_gen = lam_max(A60^{-1/2} L A60^{-1/2})
A60_inv_half = evec60 @ np.diag(1.0 / np.sqrt(ev60)) @ evec60.T
mu_gen = float(np.linalg.eigvalsh(
    A60_inv_half @ loss_mat @ A60_inv_half)[-1])

OUT['validation_vs_run12'] = {
    'lambda_min_M11': float(np.linalg.eigvalsh(M[:12, :12])[0]),
    'lambda_min_M11_run12': R12['lambda_min_M11_mine'],
    'SB60': SB60, 'SB60_run12': R12['sweep_bordered']['60']
        ['lambda_min_schur_bordered'],
    'SB800': float(np.linalg.eigvalsh(
        St[:12, :12] - St[:12, 12:] @ np.linalg.solve(
            St[12:, 12:], St[12:, :12]))[0]),
    'lambda_max_L': lam_max_L,
    'lambda_max_L_run12': R12['loss_72_to_812']['op_norm'],
    'lambda_min_SZeff': lam_min_A,
    'lambda_min_SZeff_run12': R12['floors']['lambda_min_SZeff_73_812'],
    'Lambda60_via_L': full_effect,
    'Lambda60_run12': R12['Lambda60_bordered_measured']}
print('validation vs run 12:', OUT['validation_vs_run12'])

# ---------------------------------------------------------- thresholds
# Derived BEFORE the floors: what c_floor must exceed.
colnorm = np.linalg.norm(v_tilde, axis=0)
omZ = om_all[72:]
SigmaF = float(np.sum(colnorm ** 2))
y_v0 = v_tilde.T @ v0                                # reached vector
SigmaV0 = float(y_v0 @ y_v0)
Cv = float(np.max(colnorm * omZ ** 3))
sum_om6_beyond_812 = float((4 * np.pi) ** -6
                           * (1.0 / (5 * 812.0 ** 5) + 1.0 / 812.0 ** 6))
OUT['thresholds'] = {
    'beta0': BETA0, 'beta1': BETA1,
    'SigmaF_sum_colnorm2_73_812': SigmaF,
    'SigmaV0_norm_vtildeT_v0_sq': SigmaV0,
    'c_star_Frob': SigmaF / BETA0,
    'c_star_v0': SigmaV0 / BETA0,
    'c_star_v0_vs_beta1': SigmaV0 / BETA1,
    'Cv_envelope': Cv,
    'sum_omega^-6_beyond_812': sum_om6_beyond_812,
    'c_star_tail_filed_route': 12 * Cv ** 2 * sum_om6_beyond_812 / BETA0,
    'lambda_max_L': lam_max_L,
    'v0Lv0': v0Lv0,
    'gain_measured_lammaxL_over_v0Lv0': lam_max_L / v0Lv0}
print('THRESHOLDS:', OUT['thresholds'])

# ---------------------------------------------------------- floors
def floors_at(J):
    """Floors with Z truncated at mode J (Z_J = modes 73..J)."""
    n = J - 72
    vt = v_tilde[:, :n]
    A = SZ_eff[:n, :n]
    A = 0.5 * (A + A.T)
    lamA = float(np.linalg.eigvalsh(A)[0])
    _, sv, Vt = np.linalg.svd(vt, full_matrices=False)
    U = Vt.T[:, sv > 1e-14 * sv[0]]       # Z-space basis of range(vt^T)
    comp = U.T @ A @ U
    lam_comp = float(np.linalg.eigvalsh(0.5 * (comp + comp.T))[0])
    Ainv_U = np.linalg.solve(A, U)
    inv_form = U.T @ Ainv_U
    inv_form = 0.5 * (inv_form + inv_form.T)
    c_inv = float(1.0 / np.linalg.eigvalsh(inv_form)[-1])
    LJ = vt @ np.linalg.solve(A, vt.T)
    LJ = 0.5 * (LJ + LJ.T)
    v0LJv0 = float(v0 @ LJ @ v0)
    return {'J': J, 'rank_reached': int(U.shape[1]),
            'singvals_vtilde_min': float(sv[-1]),
            'singvals_vtilde_max': float(sv[0]),
            'lambda_min_A': lamA,
            'lambda_min_compression_UTAU': lam_comp,
            'c_inv_schur_floor': c_inv,
            'v0_LJ_v0': v0LJv0,
            'c_dir_SigmaV0_over_v0LJv0': SigmaV0 / v0LJv0,
            'SigmaF_J': float(np.sum(np.linalg.norm(vt, axis=0) ** 2)),
            'bound_Frob_over_cinv': SigmaF / c_inv,
            'bound_v0_over_cinv': SigmaV0 / c_inv,
            'room_Frob_vs_beta0': BETA0 / (SigmaF / c_inv),
            'room_v0_vs_beta0': BETA0 / (SigmaV0 / c_inv)}

OUT['floors_by_truncation'] = {}
for J in (412, 612, 812):
    f = floors_at(J)
    OUT['floors_by_truncation'][str(J)] = f
    print('floors at J =', J, ':', f)

OUT['stage60'] = {
    'SB60': SB60, 'gap_lambda2_minus_lambda1': gap60,
    'lambda_max_L': lam_max_L, 'v0Lv0': v0Lv0,
    'Lambda60_full_effect': full_effect,
    'residual_norm_Lv0_minus_v0Lv0_v0': resid_Lv0,
    'mu_gen': mu_gen, 'mu_gen_times_SB60': mu_gen * SB60,
    'v0': [float(x) for x in v0]}
print('stage60:', {k: v for k, v in OUT['stage60'].items() if k != 'v0'})

# --------------------------------------- near-null directions of A
evA, evecA = np.linalg.eigh(SZ_eff)
bands = [(73, 172), (173, 372), (373, 672), (673, 812)]
def band_mass(z, a, b):
    return float(np.sum(z[a - 73:b - 72] ** 2))
_, sv812, Vt812 = np.linalg.svd(v_tilde, full_matrices=False)
U812 = Vt812.T[:, sv812 > 1e-14 * sv812[0]]
nearnull = {}
for k in range(8):
    z = evecA[:, k]
    nearnull[f'eig{k}'] = {
        'eigenvalue': float(evA[k]),
        'band_masses': {f'{a}-{b}': band_mass(z, a, b) for a, b in bands},
        'argmax_mode': int(np.argmax(np.abs(z))) + 73,
        'overlap2_with_range_vtildeT': float(np.sum((U812.T @ z) ** 2)),
        'overlap2_with_y_v0': float((z @ (y_v0 / np.sqrt(SigmaV0))) ** 2)}
OUT['nearnull_A_bottom8'] = nearnull
print('near-null bottom-8:', nearnull)
OUT['D812_diag_excess_filed'] = R12['diag_excess']['D812']
OUT['A_diag_last8'] = [float(SZ_eff[-8 + i, -8 + i]) for i in range(8)]
OUT['St_diag_last8'] = [float(St[Z, Z][-8 + i, -8 + i]) for i in range(8)]

# --------------------------------------- reached profile of y_v0, U
yhat = y_v0 / np.sqrt(SigmaV0)
OUT['reached_profile'] = {
    'y_v0_band_masses': {f'{a}-{b}': band_mass(yhat, a, b)
                         for a, b in bands},
    'y_v0_mass_beyond_412': float(np.sum(yhat[340:] ** 2)),
    'y_v0_envelope_max_abs_yj_times_omj3': float(
        np.max(np.abs(y_v0) * omZ ** 3)),
    'U_band_mass_trace_frac': {
        f'{a}-{b}': float(np.sum((U812[a - 73:b - 72, :] ** 2)))
        / U812.shape[1] for a, b in bands},
    'principal_angle_cos_max_bottom8': float(np.max(np.abs(
        U812.T @ evecA[:, :8])))}
print('reached profile:', OUT['reached_profile'])

# --------------------------------------- proof-side quantities
diagA = np.diag(SZ_eff)
prof = diagA * omZ ** 4 / np.log(omZ)
Dnorm = SZ_eff / np.sqrt(np.outer(diagA, diagA))
ev_norm = np.linalg.eigvalsh(0.5 * (Dnorm + Dnorm.T))
rho_Z = float(np.max(np.abs(ev_norm - 1.0)))
SigmaD = float(np.sum(colnorm ** 2 / diagA))
gram_vv = v_tilde @ v_tilde.T
lam_min_gram_raw = float(np.linalg.eigvalsh(0.5 * (gram_vv + gram_vv.T))[0])
c_inv_812 = OUT['floors_by_truncation']['812']['c_inv_schur_floor']
OUT['proof_quantities'] = {
    'A_diag_profile_om4_over_logom_min': float(np.min(prof)),
    'A_diag_profile_om4_over_logom_max': float(np.max(prof)),
    'A_diag_at_73': float(diagA[0]), 'A_diag_at_412': float(diagA[340]),
    'A_diag_at_812': float(diagA[-1]),
    'St_diag_at_73': float(St[Z, Z][0, 0]),
    'Ycorr_diag_rel_max': float(np.max(
        (np.diag(St[Z, Z]) - diagA) / np.diag(St[Z, Z]))),
    'Ycorr_diag_rel_at_73': float(
        (St[Z, Z][0, 0] - diagA[0]) / St[Z, Z][0, 0]),
    'Ycorr_diag_rel_at_812': float(
        (St[Z, Z][-1, -1] - diagA[-1]) / St[Z, Z][-1, -1]),
    'rho_Z_Dnormalised_offdiag_opnorm': rho_Z,
    'SigmaD_sum_colnorm2_over_Ajj': SigmaD,
    'chain_via_gram_status': 'VOID: sigma_min(v~)^2 = %.3e is below the '
        'float64 resolution (~1e-35) of the Gram matrix v~ v~^T '
        '(raw eigvalsh %.3e); the chain through lambda_min(v~ v~^T) '
        'is killed by conditioning and replaced by the profile route '
        'of Part G' % (float(sv812[-1] ** 2), lam_min_gram_raw),
    'sigma_min_vtilde': float(sv812[-1]),
    'c_inv_measured': c_inv_812}
print('proof quantities:', OUT['proof_quantities'])

# ============================================================ Part G
# The profile route: the consumed direction y' = v~^T v0 through the
# D-normalised form of A.  A = D^{1/2}(I + R)D^{1/2} with
# ||R|| = rho_Z < 1 gives the sandwich
#   (1/(1+rho_Z)) sum_j y_j^2/A_jj <= <y, A^{-1} y>
#       <= (1/(1-rho_Z)) sum_j y_j^2/A_jj .
RQ_D = float(np.sum(y_v0 ** 2 / diagA)) / SigmaV0     # per unit ||y'||^2
sand_lo = RQ_D * SigmaV0 / (1 + rho_Z)
sand_hi = RQ_D * SigmaV0 / (1 - rho_Z)
OUT['profile_route'] = {
    'RQ_D_sum_y2_over_Ajj_per_unit_norm': RQ_D,
    'c_dir_lower_bound_via_profile': (1 - rho_Z) / RQ_D,
    'v0Lv0_sandwich_lo': sand_lo, 'v0Lv0_sandwich_hi': sand_hi,
    'v0Lv0_measured': v0Lv0,
    'bound_v0Lv0_profile_route': sand_hi,
    'room_profile_route_vs_beta0': BETA0 / sand_hi,
    'room_profile_route_vs_beta1': BETA1 / sand_hi}
print('profile route:', OUT['profile_route'])

# eigen-band split of v0Lv0: where the loss actually comes from
coeff = evecA.T @ y_v0
contrib = coeff ** 2 / evA
band_edges = [(0.0, 1e-14), (1e-14, 1e-12), (1e-12, 1e-10),
              (1e-10, np.inf)]
OUT['v0Lv0_eigenband_split'] = {
    f'{lo:.0e}-{hi:.0e}': float(np.sum(contrib[(evA >= lo) & (evA < hi)]))
    for lo, hi in band_edges}
OUT['v0Lv0_eigenband_split']['total'] = float(np.sum(contrib))
OUT['v0Lv0_bottom8_contrib'] = float(np.sum(contrib[:8]))
print('v0Lv0 eigenband split:', OUT['v0Lv0_eigenband_split'])

# the direction responsible for the small UNIFORM floor: bottom
# eigenvector of the compression U^T A U, lifted back to Z-space.
comp812 = U812.T @ SZ_eff @ U812
evc, evec_c = np.linalg.eigh(0.5 * (comp812 + comp812.T))
z_worst = U812 @ evec_c[:, 0]
OUT['worst_compression_direction'] = {
    'rayleigh_quotient': float(evc[0]),
    'band_masses': {f'{a}-{b}': band_mass(z_worst, a, b)
                    for a, b in bands},
    'argmax_mode': int(np.argmax(np.abs(z_worst))) + 73,
    'overlap2_with_y_v0': float((z_worst @ yhat) ** 2)}
print('worst compression direction:', OUT['worst_compression_direction'])

# --------------------------------------- drift check (64 nodes)
print('assembling M at 64 nodes/period (drift check) ...')
M2, _ = assemble_M(64)
St2 = M2 - np.outer(u_all, u_all) / G_VAL
SY2 = St2[Y, Y]
vt2 = St2[X, Z] - St2[X, Y] @ np.linalg.solve(SY2, St2[Y, Z])
A2 = St2[Z, Z] - St2[Z, Y] @ np.linalg.solve(SY2, St2[Y, Z])
A2 = 0.5 * (A2 + A2.T)
_, sv2, Vt2 = np.linalg.svd(vt2, full_matrices=False)
U2 = Vt2.T[:, sv2 > 1e-14 * sv2[0]]
comp2 = float(np.linalg.eigvalsh(
    0.5 * (U2.T @ A2 @ U2 + (U2.T @ A2 @ U2).T))[0])
invf2 = U2.T @ np.linalg.solve(A2, U2)
c_inv2 = float(1.0 / np.linalg.eigvalsh(0.5 * (invf2 + invf2.T))[-1])
A60_2 = St2[X, X] - St2[X, Y] @ np.linalg.solve(SY2, St2[Y, X])
A60_2 = 0.5 * (A60_2 + A60_2.T)
ev60_2, evec60_2 = np.linalg.eigh(A60_2)
v0_2 = evec60_2[:, 0]
L2 = vt2 @ np.linalg.solve(A2, vt2.T)
L2 = 0.5 * (L2 + L2.T)
y_v0_2 = vt2.T @ v0_2
OUT['res64_drift'] = {
    'lambda_min_A': float(np.linalg.eigvalsh(A2)[0]),
    'lambda_min_compression': comp2,
    'c_inv': c_inv2,
    'c_inv_drift_abs': float(abs(c_inv2 - c_inv_812)),
    'compression_drift_abs': float(abs(
        comp2 - OUT['floors_by_truncation']['812']
        ['lambda_min_compression_UTAU'])),
    'c_dir_res64': float((y_v0_2 @ y_v0_2) / (v0_2 @ L2 @ v0_2)),
    'c_dir_48nodes': OUT['floors_by_truncation']['812']
        ['c_dir_SigmaV0_over_v0LJv0'],
    'v0Lv0_res64': float(v0_2 @ L2 @ v0_2),
    'SigmaV0_res64': float(y_v0_2 @ y_v0_2)}
print('res64 drift:', OUT['res64_drift'])

OUT['runtime_seconds'] = time.time() - T_START
with open('/home/hatch/workspace/gppdiscovery3/research/'
          'flint_reached_floor_2026_11_15_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('done', OUT['runtime_seconds'])
