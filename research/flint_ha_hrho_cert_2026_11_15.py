#!/usr/bin/env python
"""Flint run 16 (2026-11-15): interval certification of (H-A) and (H-rho),
the two remaining finite hypotheses for proof-grade OddFace(0.5), as
catalogued in runs 13/14 (flint_reached_floor_2026_11_15 and
flint_reached_floor_tail_2026_11_15).

STATEMENTS (run 13, section 2 (P5); L = 0.5, l = 1/4, omega_j = 4 pi j;
St = M - u u^T/g the bordered Gram on sine modes phi_j; X = 1..12,
Y = 13..72, Z = 73..812; A = S_Z^eff = St_ZZ - St_ZY St_YY^{-1} St_YZ):
  (H-A)  A_jj >= c_A * log(omega_j) / omega_j^4 for all j >= 73,
         catalogued c_A = 1.7372355768219303 (the float64 profile
         minimum of run 13, attained at j = 74).
  (H-rho) rho_Z = || D^{-1/2} A D^{-1/2} - I || <= 0.34,
         D = diag(A); float64 value 0.33419022381006247.

METHOD (all Arb balls; no float64 result is consumed anywhere):
  * Modes 1..72 and u, g: Heart's Round-7 Arb shards (GPPDiscovery5
    1da974d, experiments/round7_turbo), the same committed balls run 15
    consumed.
  * Q(omega) = R(omega) + Pi(omega): Heart's R12 closed form, the exact
    code path of run 15 (validated there against the t2rig tables).
    Master: K(omega_j) := Q(omega_j)/(pi omega_j) satisfies
        B(psi_j, psi_k) = 128 (-1)^{j+k} (K_j - K_k)
                          / (omega_j omega_k (omega_j^2 - omega_k^2))
                          - 2 bI2^2 sigma_psi,j sigma_psi,k   (j != k),
    with psi_j = phi_j - w_j delta_l, sigma_psi,j = sig(j) - w_j sinh(l/2).
    VALIDATED in-run against the shard balls at five off-diagonal pairs
    (formula balls strictly inside shard balls) -- section V1.
  * Diagonal: B(psi_j, psi_j) = (128/omega_j^2) K2_j - 2 bI2^2 sigma_psi,j^2
    with K2_j = (1/pi) int_0^inf sin^4(t/2) V(t) sin^2(t l)/(t^2-omega_j^2)^2 dt
             = (1/pi) int_{-omega_j}^{inf} sin^4(s/2) V(omega_j+s)
                      l^2 sinc^2(s l) / (2 omega_j + s)^2 ds,
    evaluated by Arb panel quadrature on s in [-omega_j, 3 omega_j] plus
    an explicit tail ball (section K2). VALIDATED in-run against the
    shard balls at j = 13, 50, 72 -- section V2.
  * St entries: St_jk = B(psi_j, psi_k) - rho_j rho_k / g for j, k in Z
    (rho_j := B(psi_j, delta_l) = (w_j/omega_j^2) q_delta + rho_j(delta),
    rho_j(delta) = w4rho_delta(j, Q_j)/omega_j^4, q_delta by run 15's
    kernel route); St_jk = B(psi_j, phi_k) - rho_j u_k/g for j in Z,
    k in Y, with B(psi_j, phi_k) = (w_j/omega_j^2) q_{phi_k}
    + w4rho_phi(k, j, Q_k, Q_j)/omega_j^4 and q_{phi_k} = H[12][k-1]
    from the shards (run 15's identity); St on Y from the shards.
  * A assembled by one Arb matrix solve. (H-A) finite range read off
    diag(A). (H-rho): Atilde = D^{-1/2} A D^{-1/2} in balls; for each
    threshold t on a ladder, the two matrices (1+t) I - Atilde and
    Atilde - (1-t) I are certified positive definite by a verified
    LDL^T: thin LDL^T of the midpoint matrix (Arb) plus the rigorous
    perturbation bound || |L^-1| R |L^-T| ||_inf < min D, R the
    entrywise radius matrix of the ball matrix. The factors for the
    catalogued threshold t = 0.34 (and the tightest certified one)
    are dumped to the .npz so the check can be rerun from the
    factors alone (verify mode).
  * (H-A) tail j > 812: analytic lower bound on B(psi_j, psi_j) from
    the s-integral form of K2 (section TAIL of the run note): core
    |s| <= 60 with the even-envelope parity cancellation, annulus
    60 < |s| <= 3000, outer regions by kernel majorants; then
    A_jj >= B(psi_j, psi_j) - rho_j^2/g - c_j with c_j the
    Y-elimination correction bounded through the fixed vector
    z = q_Y - q_delta u_Y/g and ||St_YY^{-1}||_F.

MODES:  python flint_ha_hrho_cert_2026_11_15.py k2       (compute K2 balls)
        python flint_ha_hrho_cert_2026_11_15.py certify  (full certification)
        python flint_ha_hrho_cert_2026_11_15.py verify   (re-check from dumps)
No zero data anywhere. Epistemic status: interval-certified finite
computation + one analytic tail bound, as labelled per item in the
results JSON.
"""

import glob
import io
import json
import math
import os
import subprocess
import sys
import tarfile
import time

from flint import acb, arb, arb_mat, ctx

ctx.prec = 200
T00 = time.time()

PI = arb.pi()
ELL = arb(1) / 4
SQ8 = arb(8).sqrt()
BI2 = 16 * (ELL).sinh() ** 2
SINH8 = (ELL / 2).sinh()
NS = [(2, 2), (3, 3), (4, 2), (5, 5), (7, 7), (8, 2), (9, 3), (11, 11)]
JMAX = 812
Z0 = 73                       # first Z mode
Y0, Y1 = 13, 72               # Y modes
HERE = os.path.dirname(os.path.abspath(__file__))
K2_BALLS = os.path.join(HERE, 'flint_ha_hrho_cert_2026_11_15_k2balls.json')
BALLS_OUT = os.path.join(HERE, 'flint_ha_hrho_cert_2026_11_15_balls.json')
LDLT_OUT = os.path.join(HERE, 'flint_ha_hrho_cert_2026_11_15_ldlt.npz')
RESULTS_OUT = os.path.join(HERE, 'flint_ha_hrho_cert_2026_11_15_results.json')


def ser(x):
    """Lossless ball serialization for the JSON dumps: midpoint and
    radius at 30 digits. (str(arb) prints the midpoint to ~3 digits
    and pads the radius to cover the rounding -- re-parsing it
    inflates radii up to ~4x. Caught in this run: stored K2 balls
    were up to 4x wider than the computed ones.)"""
    def num(b):
        s = b.str(30)
        return s.strip().lstrip('[').rstrip(']').split(' +/-')[0].strip()
    return f"{num(x.mid())} +/- {num(x.rad())}"


def om(j):
    return 4 * PI * j


def sig(j):
    w = om(j)
    return -SQ8 * (-1) ** j * w * SINH8 / (arb(1) / 4 + w * w)


def w_(j):
    return SQ8 * (-1) ** (j + 1) / om(j)


def sig_psi(j):
    return sig(j) - w_(j) * SINH8


# ---------------------------------------------------------------- Q(omega)
def E_(s, x, w):
    return (-s * x).exp() * (s * (w * x).sin() + w * (w * x).cos()) / (s * s + w * w)


def I_piece(a, b, w, M=110):
    tot = arb(0)
    for m in range(M + 1):
        s = arb(2 * m) + arb(1) / 2
        tot += 2 * (E_(s, a, w) - E_(s, b, w))
    s0 = 2 * (M + 1) + arb(1) / 2
    return tot + arb(0, 4 * (-s0 * a).exp() / s0 / (1 - (-2 * a).exp()))


def I_first(w, M=110):
    z = acb(arb(1) / 4) + acb(0, 1) * (w / 2)
    tot = z.digamma().imag
    for m in range(M + 1):
        s = arb(2 * m) + arb(1) / 2
        tot -= 2 * E_(s, arb(1) / 2, w)
    s0 = 2 * (M + 1) + arb(1) / 2
    return tot + arb(0, 4 * (-s0 / 2).exp() / s0 / (1 - arb(-1).exp()))


PIECES = [(arb(0), arb(1) / 2, arb(3) / 16), (arb(1) / 2, arb(1), arb(1) / 8),
          (arb(1), arb(3) / 2, -arb(1) / 8), (arb(3) / 2, arb(2), -arb(1) / 32),
          (arb(2), arb(5) / 2, arb(1) / 32)]


def R_of(w):
    tot = arb(0)
    for (a, b, s) in PIECES:
        tot += s * (I_first(w) if a == 0 else I_piece(a, b, w))
    return (PI / 2) * tot


CS = [arb(3) / 16, -arb(1) / 16, -arb(1) / 4, arb(3) / 32, arb(1) / 16, -arb(1) / 32]


def s_of(u):
    uf = float(u.mid())
    tot = arb(0)
    for k in range(6):
        if k / 2 < uf:
            tot += CS[k]
    return tot


def Pi_of(w):
    return PI * sum(((arb(p).log() / arb(n).sqrt()) * s_of(arb(n).log())
                     * (w * arb(n).log()).sin() for n, p in NS), arb(0))


def Q_of(w):
    return R_of(w) + Pi_of(w)


PREF = (16 / PI) * SQ8


def w4rho_delta(j, Qj):
    return PREF * (-1) ** j * Qj - (BI2 ** 2 / 8) * SINH8 * sig(j)


def w4rho_phi(i, j, Qi, Qj):
    wi, wj = om(i), om(j)
    return (128 / PI) * (-1) ** (i + j) * (wi * wj / (wj * wj - wi * wi)) \
        * (Qj / wj - Qi / wi) - (BI2 ** 2 / 8) * sig(i) * sig(j)


def K_of(j, Qj):
    """Master K(omega_j) = Q(omega_j)/(pi omega_j)."""
    return Qj / (PI * om(j))


# ---------------------------------------------------------------- shards
def load_shards():
    d = '/tmp/heart-src/experiments/round7_turbo'
    if glob.glob(d + '/shard_*.json'):
        sd = d
    else:
        out = '/tmp/flint16/shardx'
        if not glob.glob(out + '/experiments/round7_turbo/shard_*.json'):
            os.makedirs(out, exist_ok=True)
            data = subprocess.run(
                ['git', '-C', os.path.expanduser('~/workspace/GPPDiscovery5'),
                 'archive', '1da974d', 'experiments/round7_turbo'],
                capture_output=True, check=True).stdout
            tf = tarfile.open(fileobj=io.BytesIO(data))
            tf.extractall(out)
        sd = out + '/experiments/round7_turbo'
    B = {}
    for fn in sorted(glob.glob(sd + '/shard_*.json')):
        for k, v in json.load(open(fn)).items():
            a, b = k.split(',')
            B[('d' if a == 'd' else int(a), 'd' if b == 'd' else int(b))] = arb(v)
    return B


def shard_tables(B):
    def Bph(j, k):
        return B[(j - 1, k - 1)] if j <= k else B[(k - 1, j - 1)]

    def Bpd(j):
        return B[(j - 1, 'd')]

    g = B[('d', 'd')]
    # q_{phi_i} = H[12][i-1] = om_i^2 (B(delta, phi_i) - w_i g)   (run 15)
    q_phi = {i: om(i) ** 2 * (Bpd(i) - w_(i) * g) for i in range(1, 73)}
    return Bph, Bpd, g, q_phi


def Bpsi_shard(Bph, Bpd, g, j, k):
    return Bph(j, k) - w_(j) * Bpd(k) - w_(k) * Bpd(j) + w_(j) * w_(k) * g


# ---------------------------------------------------------------- q_delta
# Run 15's kernel route, verbatim in substance: q_delta = -1/2 [K''(0) - K''(1/2)],
# K''(d) = (BI2^2/2) cosh(d/2) - prime2(d) + arch2(d).
LOGN = {n: arb(n).log() for n, p in NS}


def beta2_val(y):
    a = abs(y)
    alo = float(a.lower())
    ahi = float(a.upper())

    def f1(x):
        return -2 + 3 * x

    def f2(x):
        return 2 - x

    if ahi < 1:
        return f1(a)
    if alo > 2:
        return arb(0)
    if alo > 1 and ahi < 2:
        return f2(a)
    cands = []
    if alo < 1:
        cands.append(f1(a))
    if ahi > 1 and alo < 2:
        cands.append(f2(a))
    if ahi > 2:
        cands.append(arb(0))
    lo = min(float(c.lower()) for c in cands)
    hi = max(float(c.upper()) for c in cands)
    return arb(lo) + arb(0, hi - lo) if hi > lo else cands[0]


def beta2_lin(c, sgn, umid):
    ymid = c + sgn * umid
    a = abs(float(ymid.mid()))
    if a > 2:
        return (arb(0), arb(0))
    if a < 1:
        if float((c + sgn * umid).mid()) >= 0:
            return (-2 + 3 * c, 3 * sgn)
        return (-2 - 3 * c, -3 * sgn)
    if float((c + sgn * umid).mid()) >= 0:
        return (2 - c, -sgn)
    return (2 + c, sgn)


MK = 130


def H0s(a, b):
    tot = arb(0)
    for m in range(MK + 1):
        s = arb(2 * m) + arb(1) / 2
        tot += 2 * ((-s * a).exp() - (-s * b).exp()) / s
    s0 = arb(2 * (MK + 1)) + arb(1) / 2
    return tot + arb(0, 2 * (-s0 * a).exp() / (s0 * (1 - (-2 * a).exp())))


def H1s(a, b):
    tot = arb(0)
    for m in range(MK + 1):
        s = arb(2 * m) + arb(1) / 2
        tot += 2 * ((-s * a).exp() * (s * a + 1) - (-s * b).exp() * (s * b + 1)) / (s * s)
    s0 = arb(2 * (MK + 1)) + arb(1) / 2
    return tot + arb(0, 2 * (-s0 * a).exp() * (a + 2 / s0) / (s0 * (1 - (-2 * a).exp())))


EPS = arb('1e-6')


def _integ0(f, b, lead):
    if float(b.mid()) <= 1e-6:
        return lead(b)
    return lead(EPS) + acb.integral(f, acb(EPS), acb(b)).real


def R_int(b):
    def f(u, _):
        return u * (u / 2).exp() / u.sinh() - 1 - u / 2

    def lead(x):
        return -x ** 3 / 72 - x ** 4 / 64 + arb(0, 2 * x ** 5)

    return _integ0(f, b, lead)


def R1_int(b):
    def f(u, _):
        return (u / 2).exp() / u.sinh() - 1 / u - arb(1) / 2

    def lead(x):
        return -x ** 2 / 48 - x ** 3 / 48 + arb(0, 2 * x ** 4)

    return _integ0(f, b, lead)


def H0f(a, b):
    if float(a.mid()) < 0.5:
        return (b / a).log() + (b - a) / 2 + (R1_int(b) - R1_int(a))
    return H0s(a, b)


def H1f(a, b):
    if float(a.mid()) < 0.5:
        return (b - a) + (b * b - a * a) / 4 + (R_int(b) - R_int(a))
    return H1s(a, b)


def H1_0(b):
    return b + b * b / 4 + R_int(b)


_H03 = None
_W0 = None


def arch2(d):
    global _H03, _W0
    if _H03 is None:
        _H03 = sum((2 * (-3 * (arb(2 * m) + arb(1) / 2)).exp() / (arb(2 * m) + arb(1) / 2)
                    for m in range(MK + 1)), arb(0)) + arb(0, 1e-150)
        _W0 = (acb(arb(1) / 4)).digamma().real - PI.log()
    b2d = beta2_val(arb(d))
    ad = abs(float(arb(d).mid()))
    pts = sorted(set([0.0, ad, 1 - ad, 1 + ad, 2 - ad, 2 + ad, 3.0]))
    pts = [p for p in pts if 0 <= p <= 3.0]
    tot = _W0 * b2d
    for a, b in zip(pts[:-1], pts[1:]):
        if b - a < 1e-14:
            continue
        aa, bb = arb(str(a)), arb(str(b))
        umid = (a + b) / 2
        pp, qp = beta2_lin(arb(d), 1, umid)
        pm, qm = beta2_lin(arb(d), -1, umid)
        alpha = b2d - (pp + pm) / 2
        beta = -(qp + qm) / 2
        if a == 0.0:
            tot += beta * H1_0(bb)
        else:
            tot += alpha * H0f(aa, bb) + beta * H1f(aa, bb)
    return tot + b2d * _H03


def prime2(d):
    dd = arb(d)
    return sum(((arb(p).log() / arb(n).sqrt())
                * (beta2_val(LOGN[n] + dd) + beta2_val(dd - LOGN[n]))
                for n, p in NS), arb(0))


def K2pp(d):
    """K_S''(d) in balls (run 15 naming: K2 there)."""
    dd = arb(d)
    return (BI2 ** 2 / 2) * (dd / 2).cosh() - prime2(dd) + arch2(dd)


def q_delta_ball():
    return -(K2pp(arb(0)) - K2pp(arb(1) / 2)) / 2


# ---------------------------------------------------------------- K2 diag
# K2_j = (1/pi) int_{-om}^{inf} sin^4(s/2) V(om+s) l^2 sinc^2(s l)
#          / (2 om + s)^2 ds,  V(t) = Re psi(1/4 + i t/2) - log pi
#          - 2 sum_n (Lam(n)/sqrt n) cos(t log n).
# Rigorous bound |V(t)| <= log t + 12 for t >= 1, |V| <= 23 on [0, 1]:
# Binet: |Re psi(1/4+it/2) - log(t/2)| <= 1.5 for t >= 1 (the alternating
# Binet remainder for Re psi(x+iy) with x = 1/4 is < 1/(2(x^2+y^2)) +
# int bound <= 1.5 in fact far less; we consume 1.5), |V| <= log(t/2)
# + 1.5 + log pi + 2 P0 with P0 = sum Lam/sqrt n = 4.2613..., so
# |V| <= log t + 0.807 + 8.523 < log t + 12. On [0,1]: Re psi(1/4 + it/2)
# is decreasing in |t| from psi(1/4) = -gamma - pi/2 - 3 log 2, so
# |Re psi| <= 4.228 and |V| <= 4.228 + 1.145 + 8.523 < 23.
LOGN_A = {n: arb(n).log() for n, p in NS}
LAMN = {n: arb(p).log() / arb(n).sqrt() for n, p in NS}


def Vval(t):
    z = acb(arb(1) / 4) + acb(0, 1) * (t / 2)
    v = z.digamma().real - PI.log()
    for n, p in NS:
        v -= 2 * LAMN[n] * (t * LOGN_A[n]).cos()
    return v


def sinc2_l(s):
    """l^2 sinc^2(s l) = sin^2(s l)/s^2, ball-safe at s = 0 (limit l^2)."""
    x = s * ELL
    if abs(float(x.mid())) < 0.05:
        # sin(x)/x = 1 - x^2/6 + x^4/120 - ... ; |x| < 0.05: tail < x^6/5040
        x2 = x * x
        r = 1 - x2 / 6 + x2 * x2 / 120
        r += arb(0, float((abs(x2) ** 3 / 5040).upper()) + 1e-30)
        return (ELL * r) ** 2
    return (x.sin() / s) ** 2


def K2_integrand(s_acb, w):
    s = s_acb.real
    sn = (s / 2).sin()
    return acb((sn ** 4) * Vval(w + s) * sinc2_l(s) / (2 * w + s) ** 2)


def K2_ball(w, step=32):
    """Ball for K2 at frequency w = omega_j, panels on s in [-w, 3w]
    plus an explicit tail ball for s > 3w."""
    old = ctx.prec
    ctx.prec = 160
    try:
        wm = float(w.mid())
        lo = -wm
        hi = 3 * wm
        acc = arb(0)
        nfull = int((hi - lo) / step)
        for i in range(nfull):
            a = lo + i * step
            acc += acb.integral(lambda u, _: K2_integrand(u, w),
                                acb(a), acb(a + step), abs_tol=1e-15).real
        rem = hi - (lo + nfull * step)
        if rem > 0:
            acc += acb.integral(lambda u, _: K2_integrand(u, w),
                                acb(lo + nfull * step), acb(hi),
                                abs_tol=1e-15).real
        # tail s > 3w: |integrand| <= (log(w+s)+12)/(s^2 (2w+s)^2);
        # split at s = 6w: on [3w, 6w], log(w+s) <= log(7w) and
        # (2w+s)^2 >= (5w)^2; on [6w, inf), (2w+s)^2 >= s^2 and
        # w + s <= 2s.
        lw = arb(wm).log()
        t1 = (lw + arb(7).log() + 12) / (75 * w ** 3)
        t2 = (((12 * w).log() + 12) / 3 + arb(1) / 9) / (6 * w) ** 3
        tail = (t1 + t2) / PI
        return acc / PI + arb(0, float(tail.upper()))
    finally:
        ctx.prec = old


def Bpsi_diag_formula(j, K2j):
    return (128 / om(j) ** 2) * K2j - 2 * BI2 ** 2 * sig_psi(j) ** 2


def Bpsi_offdiag(j, k, Kj, Kk):
    a = om(j) ** 2 - om(k) ** 2
    return 128 * (-1) ** (j + k) * (Kj - Kk) / (om(j) * om(k) * a) \
        - 2 * BI2 ** 2 * sig_psi(j) * sig_psi(k)


def validate_shards(Bph, Bpd, g, Q):
    """V1: off-diagonal formula vs shards (5 pairs); V2: K2 diag vs shards."""
    rep = {}
    ok = True
    for (j, k) in [(13, 14), (13, 72), (20, 50), (50, 72), (31, 47)]:
        fb = Bpsi_offdiag(j, k, K_of(j, Q[j]), K_of(k, Q[k]))
        sb = Bpsi_shard(Bph, Bpd, g, j, k)
        inside = float(fb.lower()) >= float(sb.lower()) - 1e-30 and \
            float(fb.upper()) <= float(sb.upper()) + 1e-30
        rep[f'V1_{j}_{k}'] = {'formula': str(fb), 'shard': str(sb),
                              'formula_inside_shard': inside,
                              'formula_rad': float(fb.rad())}
        ok = ok and inside
    for j in (13, 50, 72):
        K2j = K2_ball(om(j))
        fb = Bpsi_diag_formula(j, K2j)
        sb = Bpsi_shard(Bph, Bpd, g, j, j)
        ov = fb.overlaps(sb)
        rep[f'V2_{j}'] = {'formula': str(fb), 'shard': str(sb),
                          'K2': str(K2j), 'overlap': bool(ov),
                          'formula_relrad': float((fb.rad() / abs(fb)).upper())
                          if fb != 0 else None}
        ok = ok and bool(ov)
    rep['V_all_ok'] = ok
    return rep


def compute_k2_all():
    """K2 balls for j = Z0..JMAX, checkpointed to K2_BALLS (resumable)."""
    done = {}
    if os.path.exists(K2_BALLS):
        done = {int(k): arb(v) for k, v in json.load(open(K2_BALLS)).items()}
    todo = [j for j in range(Z0, JMAX + 1) if j not in done]
    print(f'K2: {len(done)} done, {len(todo)} to go', flush=True)
    for idx, j in enumerate(todo):
        done[j] = K2_ball(om(j))
        if idx % 20 == 19 or idx == len(todo) - 1:
            json.dump({str(k): ser(v) for k, v in sorted(done.items())},
                      open(K2_BALLS, 'w'))
            print(f'K2: {len(done)}/{JMAX - Z0 + 1} '
                  f'({time.time() - T00:.0f}s)', flush=True)
    json.dump({str(k): ser(v) for k, v in sorted(done.items())},
              open(K2_BALLS, 'w'))
    return done


# ---------------------------------------------------------------- assembly
def rho_delta_ball(j, Qj):
    """rho_j(delta_l) = B(psi_j, delta) - (w_j/om_j^2) q_delta part removed:
    precisely rho_j := B(psi_j, delta_l) = u_j - w_j g and the closed form
    rho_j(delta) of run 13 is w4rho_delta(j, Q_j)/om_j^4; the full rho_j
    used in St is (w_j/om_j^2) q_delta + rho_j(delta) (run 13, (P3))."""
    return w4rho_delta(j, Qj) / om(j) ** 4


def assemble(Bph, Bpd, g, q_phi, qd, Q, K2):
    """Returns dict with A (arb_mat 740x740), StYY, and scalar balls."""
    Ys = list(range(Y0, Y1 + 1))
    Zs = list(range(Z0, JMAX + 1))
    ny, nz = len(Ys), len(Zs)
    u = {j: Bpd(j) for j in range(1, 73)}
    # St on Y (shards)
    StYY = arb_mat([[Bph(a, b) - u[a] * u[b] / g for b in Ys] for a in Ys])
    # rho_j for j in Z
    rho = {j: w_(j) / om(j) ** 2 * qd + rho_delta_ball(j, Q[j]) for j in Zs}
    K = {j: K_of(j, Q[j]) for j in range(1, JMAX + 1)}
    # St_ZZ
    rows = []
    for j in Zs:
        rj = rho[j]
        row = []
        for k in Zs:
            if j == k:
                b = Bpsi_diag_formula(j, K2[j])
            else:
                b = Bpsi_offdiag(j, k, K[j], K[k])
            row.append(b - rj * rho[k] / g)
        rows.append(row)
    StZZ = arb_mat(rows)
    # St_ZY: B(psi_j, phi_k) - rho_j u_k / g
    rows = []
    for j in Zs:
        rj = rho[j]
        wj = w_(j)
        row = []
        for k in Ys:
            b = wj / om(j) ** 2 * q_phi[k] \
                + w4rho_phi(k, j, Q[k], Q[j]) / om(j) ** 4
            row.append(b - rj * u[k] / g)
        rows.append(row)
    StZY = arb_mat(rows)
    StYZ = StZY.transpose()
    X = StYY.solve(StYZ)                    # 60 x 740
    A = StZZ - StZY * X
    # scalars for the tail bound: C_z = z^T StYY^{-1} z, z = q_Y - qd u_Y/g
    zvec = arb_mat([[q_phi[k] - qd * u[k] / g] for k in Ys])
    Cz = (zvec.transpose() * StYY.solve(zvec))[0, 0]
    invYY = StYY.inv()
    fro2 = arb(0)
    for i in range(ny):
        for jj in range(ny):
            e = invYY[i, jj]
            fro2 += (abs(e.mid()) + e.rad()) ** 2
    invF_up = float(fro2.sqrt().upper())
    uY2 = arb(0)
    for k in Ys:
        e = u[k]
        uY2 += (abs(e.mid()) + e.rad()) ** 2
    return {'A': A, 'StYY': StYY, 'Cz': Cz, 'invF_up': invF_up,
            'uY_up': float(uY2.sqrt().upper()), 'rho': rho, 'X': X}


# ---------------------------------------------------------------- LDL^T
def ldlt_full(M, block=92):
    """Block LDL^T of a symmetric arb_mat. Returns (L, D) with
    M = L diag(D) L^T (in balls), or None if a pivot fails."""
    n = M.nrows()
    Lfull = [[arb(0)] * n for _ in range(n)]
    Dfull = [arb(0)] * n
    W = M
    off = 0
    while True:
        m = W.nrows()
        e = min(block, m)
        Db = [[W[i, j] for j in range(e)] for i in range(e)]
        fac = _scalar_ldlt(Db)
        if fac is None:
            return None
        D, L = fac
        for i in range(e):
            Dfull[off + i] = D[i]
            for j in range(e):
                Lfull[off + i][off + j] = L[i][j]
        if e == m:
            break
        Lm = arb_mat(L)
        A21 = arb_mat([[W[i, j] for j in range(e)] for i in range(e, m)])
        # Y = A21 * L^{-T}:  Y^T = solve(L, A21^T)
        Y = Lm.solve(A21.transpose()).transpose()
        Ys = arb_mat([[Y[i, j] / D[j] for j in range(e)]
                      for i in range(m - e)])
        for i in range(m - e):
            for j in range(e):
                Lfull[off + e + i][off + j] = Ys[i, j]
        W22 = arb_mat([[W[i, j] for j in range(e, m)] for i in range(e, m)])
        W = W22 - Ys * Y.transpose()
        off += e
    return arb_mat(Lfull), Dfull


def certify_pd(P, R):
    """Verified PD certificate for a ball matrix P with entrywise radius
    matrix R (floats, upper bounds): thin LDL^T of the midpoint matrix
    plus the perturbation bound || |L^-1| R |L^-T| ||_inf < min D.
    Returns dict with the certificate data, or {'ok': False}."""
    import numpy as np
    n = P.nrows()
    Mthin = arb_mat([[P[i, k].mid() for k in range(n)] for i in range(n)])
    out = ldlt_full(Mthin)
    if out is None:
        return {'ok': False, 'stage': 'thin_ldlt'}
    L, D = out
    dmin = min(float(x.lower()) for x in D)
    Linv = L.inv()
    Ainv = np.array([[float((abs(Linv[i, k].mid()) + Linv[i, k].rad()).upper())
                      for k in range(n)] for i in range(n)])
    Wm = Ainv @ R @ Ainv.T
    rowmax = float(Wm.max(axis=1).max()) * (1 + 1e-6)
    return {'ok': bool(rowmax < dmin), 'dmin': dmin, 'W_rowmax': rowmax,
            'L': L, 'D': D}


def check_ldlt_dump(At, R_at):
    """Independent check of the dumped LDL^T factors (verify mode):
    rebuilds P for the dumped thresholds, forms S = M - L D L^T in
    Arb, and re-runs the perturbation bound with R_total = R + |S|."""
    import numpy as np
    nz = At.nrows()
    npz = np.load(LDLT_OUT)
    out = {}
    for key in ('cat', 'tight'):
        if f'{key}_L1_hi' not in npz:
            continue
        t = float(npz[f'{key}_t'][0])
        res = {}
        for sfx, sgn in (('1', -1), ('2', +1)):
            Lhi = npz[f'{key}_L{sfx}_hi']
            Llo = npz[f'{key}_L{sfx}_lo']
            Dhi = npz[f'{key}_D{sfx}_hi']
            Dlo = npz[f'{key}_D{sfx}_lo']
            L = arb_mat([[arb(float(Lhi[i, k])) + arb(float(Llo[i, k]))
                          for k in range(nz)] for i in range(nz)])
            Dv = [arb(float(Dhi[i])) + arb(float(Dlo[i]))
                  for i in range(nz)]
            if sgn < 0:
                P = arb_mat([[((arb(1) + arb(str(t)) - At[i, i]) if i == k
                               else -At[i, k]) for k in range(nz)]
                             for i in range(nz)])
            else:
                P = arb_mat([[((At[i, i] - (arb(1) - arb(str(t)))) if i == k
                               else At[i, k]) for k in range(nz)]
                             for i in range(nz)])
            Mthin = arb_mat([[P[i, k].mid() for k in range(nz)]
                             for i in range(nz)])
            Dm = arb_mat([[Dv[i] if i == j else arb(0) for j in range(nz)]
                          for i in range(nz)])
            S = Mthin - L * Dm * L.transpose()
            Sabs = np.array([[float((abs(S[i, k].mid()) + S[i, k].rad())
                                     .upper()) for k in range(nz)]
                             for i in range(nz)])
            Linv = L.inv()
            Ainv = np.array([[float((abs(Linv[i, k].mid())
                                     + Linv[i, k].rad()).upper())
                              for k in range(nz)] for i in range(nz)])
            Wm = Ainv @ (R_at + Sabs) @ Ainv.T
            rowmax = float(Wm.max(axis=1).max()) * (1 + 1e-6)
            dmin = min(float(Dhi[i]) + float(Dlo[i]) for i in range(nz))
            res[sfx] = {'dmin': dmin, 'W_rowmax': rowmax,
                        'S_max': float(Sabs.max()), 'ok': bool(rowmax < dmin)}
        out[key] = {'t': t, 'ok': bool(res['1']['ok'] and res['2']['ok']),
                    'parts': res}
    return out


def _scalar_ldlt(M):
    """M: list of lists of arb (symmetric). Returns (pivots, L, D) or None."""
    n = len(M)
    L = [[arb(0)] * n for _ in range(n)]
    D = [arb(0)] * n
    for i in range(n):
        L[i][i] = arb(1)
    for j in range(n):
        s = M[j][j]
        for k in range(j):
            s -= L[j][k] ** 2 * D[k]
        if not (float(s.lower()) > 0):
            return None
        D[j] = s
        for i in range(j + 1, n):
            t = M[i][j]
            for k in range(j):
                t -= L[i][k] * L[j][k] * D[k]
            L[i][j] = t / s
    return D, L


def block_ldlt_min_pivot(M, block=92):
    """Block LDL^T of a symmetric arb_mat. Returns (ok, min pivot lower)."""
    n = M.nrows()
    W = M
    minpiv = None
    pivots = []
    s = 0
    while s < n:
        e = min(s + block, n)
        Db = [[W[i, j] for j in range(s, e)] for i in range(s, e)]
        fac = _scalar_ldlt(Db)
        if fac is None:
            return False, None, pivots
        D, L = fac
        for p in D:
            lp = float(p.lower())
            pivots.append(lp)
            minpiv = lp if minpiv is None else min(minpiv, lp)
        if e < n:
            nb = e - s
            Lm = arb_mat(L)
            Dm = arb_mat([[D[i] if i == j else arb(0) for j in range(nb)]
                          for i in range(nb)])
            A21 = arb_mat([[W[i, j] for j in range(s, e)]
                           for i in range(e, n)])
            # Y = A21 * L^{-T}:  Y^T = solve(L, A21^T)
            Yt = Lm.solve(A21.transpose())
            Y = Yt.transpose()
            Ys = arb_mat([[Y[i, j] / D[j] for j in range(nb)]
                          for i in range(n - e)])
            W22 = arb_mat([[W[i, j] for j in range(e, n)]
                           for i in range(e, n)])
            W = W22 - Ys * Y.transpose()
            # rebuild full W for the remaining loop: only trailing block
            # is needed; emulate by keeping W as the trailing matrix and
            # shifting s to 0 on it.
            n = W.nrows()
            s = 0
            # NOTE: pivots already recorded for processed blocks; the
            # trailing matrix is the exact Schur complement.
            continue
        break
    return True, minpiv, pivots


# ---------------------------------------------------------------- tail
C_DELTA = arb('19.5771')          # Heart R12 (T2-rig), banked
C_PHI = arb('69.5492')            # Heart R12 (T2-rig), banked
D_ENV = [(arb(1), arb(1) / 2), (arb(4), arb(1)), (-arb(3) / 2, arb(3) / 2),
         (-arb(1), arb(2)), (arb(1) / 2, arb(5) / 2)]


def tailmass_up(S):
    """Upper bound for int_{|s|>S} E(s) ds, E = sin^4(s/2) sin^2(s/4)/s^2."""
    S = arb(S)
    return (arb(3) / S + arb('15.4') / S ** 2) / 8


def ihat(mu):
    """int_{-inf}^{inf} E(s) cos(s mu) ds = (pi/16) sum d_k (nu_k - mu)_+."""
    tot = arb(0)
    for d, nu in D_ENV:
        if float(nu.mid()) > float(mu.mid()):
            tot += d * (nu - mu)
    return (PI / 16) * tot


def wh_ball(t):
    """wh(t) = Re psi(1/4 + i t/2) - log pi, Arb ball. wh is increasing
    in t > 0: d/dt Re psi(1/4 + it/2) = -(1/2) Im psi_1(1/4 + it/2) and
    Im psi_1(x + iy) = -2xy sum_n 1/((n+x)^2 + y^2)^2 < 0 for y > 0."""
    z = acb(arb(1) / 4) + acb(0, 1) * (arb(t) / 2)
    return z.digamma().real - PI.log()


def tail_profile_lower(g_lo, qd_up, Cz_up, invF_up, uY_up):
    """Certified lower bound for A_jj * om_j^4 / log om_j, all j >= 813,
    evaluated at om = om_813 with monotone majorants (run note, TAIL)."""
    omb = float(om(813).mid())
    omB = arb(omb)
    logom = omB.log()
    P0 = sum((LAMN[n] for n, p in NS), arb(0))
    I0 = 3 * PI / 32
    tm60 = tailmass_up(60)
    core_mass_lo = I0 - tm60
    # core: wh part (wh increasing; evaluated at om_813 thresholds)
    Wlo = wh_ball(omb - 60).lower()
    Bs = 2 * arb(60).log() + arb('0.01')          # int_core E|s| ds bound
    wh_up = (omB / 2).log() + arb('0.356') + arb(120) / omB  # wh(om+60) up
    slope_corr = (wh_up * Bs) / omB
    bracket_wh = Wlo * core_mass_lo - slope_corr
    # core: prime part
    prime = arb(0)
    for n, p in NS:
        prime += LAMN[n] * (abs(ihat(LOGN_A[n])) + arb('0.001746') + tm60)
    bracket = bracket_wh - 2 * prime
    core = (32 / (PI * logom)) * bracket
    # annulus 60 < |s| <= 3000
    ann_mass_up = tailmass_up(60)
    Vlo_ann_abs = 2 * P0 - wh_ball(omb - 3000).lower()
    ann = (32 / (PI * logom)) * Vlo_ann_abs * ann_mass_up
    # right outer s > 3000
    rout = (32 / (PI * logom)) * ((3 * omB).log() + 12) * tailmass_up(3000) \
        + (16 / (PI * logom)) * (((4 * omB).log() + 12) / 3 + arb(1) / 9) / omB
    # left region t in [0, om-3000]
    fA = (omB / 2).log() + 12
    int_left1 = arb('13.9') + (arb(3) / 8) * (omB / 2) * ((omB / 2).log() + 11) \
        + (arb(1) / 2) * (fA + 12 + (omB / 2).log()) \
        + (arb(1) / 8) * (fA + 12 + (omB / 2).log())
    left1 = (128 / (PI * logom)) * (arb(16) / 9) * int_left1 / omB ** 2
    left2 = (128 / (PI * logom)) * ((omB).log() + 12) / (arb('2.25') * 3000)
    # rho_j^2 / g and Y-elimination c_j, profile units at om_813
    rho_up = SQ8 * arb(qd_up) / omB ** 3 + C_DELTA / omB ** 4
    rho_term = (rho_up ** 2 / arb(g_lo)) * omB ** 4 / logom
    sigY = arb('0.35451') / om(13)
    sigj = arb('0.35451') / omB
    rhoY_up = arb(60).sqrt() * (C_PHI / (omB - om(72))
                                + (BI2 ** 2 / 8) * sigY * sigj) / omB ** 4
    bterm = arb(invF_up) * (rhoY_up + (C_DELTA / omB ** 4)
                            * arb(uY_up) / arb(g_lo)) ** 2
    aterm = 8 * arb(Cz_up) / omB ** 6
    cj = (aterm.sqrt() + bterm.sqrt()) ** 2
    cj_term = cj * omB ** 4 / logom
    total = core - ann - rout - left1 - left2 - rho_term - cj_term
    parts = {'core': float(core.lower()), 'annulus': float(ann.upper()),
             'right_outer': float(rout.upper()), 'left1': float(left1.upper()),
             'left2': float(left2.upper()), 'rho_term': float(rho_term.upper()),
             'cj_term': float(cj_term.upper()),
             'bracket_wh': float(bracket_wh.lower()),
             'prime_bracket': float((2 * prime).upper()),
             'tail_total': float(total.lower())}
    return float(total.lower()), parts


# ---------------------------------------------------------------- driver
def ball_up(x):
    return float((abs(x.mid()) + x.rad()).upper())


def run_certification(from_dumps=False):
    B = load_shards()
    Bph, Bpd, g, q_phi = shard_tables(B)
    if from_dumps:
        dump = json.load(open(BALLS_OUT))
        Q = {int(k): arb(v) for k, v in dump['Q'].items()}
        K2 = {int(k): arb(v) for k, v in json.load(open(K2_BALLS)).items()}
    else:
        Q = {j: Q_of(om(j)) for j in range(1, JMAX + 1)}
        print(f'Q done ({time.time() - T00:.0f}s)', flush=True)
        if not os.path.exists(K2_BALLS):
            raise SystemExit('K2 balls missing: run mode "k2" first')
        K2 = {int(k): arb(v) for k, v in json.load(open(K2_BALLS)).items()}
    val = validate_shards(Bph, Bpd, g, Q)
    print('validations:', json.dumps(val, default=str)[:400], flush=True)
    if not val['V_all_ok']:
        raise SystemExit('VALIDATION FAILED -- not certifying')
    qd = q_delta_ball()
    print(f'q_delta = {qd} ({time.time() - T00:.0f}s)', flush=True)
    asm = assemble(Bph, Bpd, g, q_phi, qd, Q, K2)
    A = asm['A']
    nz = JMAX - Z0 + 1
    print(f'assembled ({time.time() - T00:.0f}s)', flush=True)
    # (H-A) finite range
    prof_lo = {}
    for idx, j in enumerate(range(Z0, JMAX + 1)):
        p = A[idx, idx] * om(j) ** 4 / om(j).log()
        prof_lo[j] = float(p.lower())
    jmin = min(prof_lo, key=prof_lo.get)
    cA_fin = prof_lo[jmin]
    # float cross-check of midpoints against run 13
    import numpy as np
    Amid = np.array([[float(A[i, k].mid()) for k in range(nz)]
                     for i in range(nz)])
    d = np.sqrt(np.diag(Amid))
    At = Amid / np.outer(d, d)
    ev = np.linalg.eigvalsh(At)
    rho_float = float(max(abs(ev - 1)))
    prof_float_min = float(np.min(np.diag(Amid)
                                  * (4 * np.pi * np.arange(Z0, JMAX + 1)) ** 4
                                  / np.log(4 * np.pi * np.arange(Z0, JMAX + 1))))
    # (H-rho): Atilde balls; verified LDL^T certificates on a ladder of
    # thresholds. rho_Z <= t  iff  (1+t) I - Atilde >= 0 and
    # Atilde - (1-t) I >= 0.
    import numpy as np
    sq = [A[i, i].sqrt() for i in range(nz)]
    At = arb_mat([[A[i, k] / (sq[i] * sq[k]) for k in range(nz)]
                  for i in range(nz)])
    R_at = np.array([[float(At[i, k].rad().upper()) for k in range(nz)]
                     for i in range(nz)])
    ladder = [0.34, 0.25, 0.20, 0.18, 0.17, 0.165, 0.16, 0.158]
    rho_runs = {}
    ldlt_arrays = {}
    for t in ladder:
        c1 = arb(1) + arb(str(t))
        c0 = arb(1) - arb(str(t))
        P1 = arb_mat([[((c1 - At[i, i]) if i == k else -At[i, k])
                       for k in range(nz)] for i in range(nz)])
        P2 = arb_mat([[((At[i, i] - c0) if i == k else At[i, k])
                       for k in range(nz)] for i in range(nz)])
        r1 = certify_pd(P1, R_at)
        r2 = certify_pd(P2, R_at)
        ok = bool(r1.get('ok') and r2.get('ok'))
        rho_runs[str(t)] = {
            'certified': ok,
            'P1_dmin': r1.get('dmin'), 'P1_W_rowmax': r1.get('W_rowmax'),
            'P2_dmin': r2.get('dmin'), 'P2_W_rowmax': r2.get('W_rowmax')}
        print(f'rho <= {t}: certified={ok} '
              f'(dmin {r1.get("dmin")}, {r2.get("dmin")}; '
              f'W {r1.get("W_rowmax")}, {r2.get("W_rowmax")}) '
              f'({time.time() - T00:.0f}s)', flush=True)
        if ok and r1.get('L') is not None and r2.get('L') is not None:
            for key in (('cat', 'tight') if t == 0.34 else ('tight',)):
                ldlt_arrays[f'{key}_t'] = np.array([t])
                for sfx, r in (('1', r1), ('2', r2)):
                    hi = np.array([[float(r['L'][i, k].mid())
                                    for k in range(nz)] for i in range(nz)])
                    lo = np.array([[float((r['L'][i, k].mid()
                                           - arb(float(r['L'][i, k].mid())))
                                          .mid()) for k in range(nz)]
                                   for i in range(nz)])
                    ldlt_arrays[f'{key}_L{sfx}_hi'] = hi
                    ldlt_arrays[f'{key}_L{sfx}_lo'] = lo
                    dhi = np.array([float(x.mid()) for x in r['D']])
                    dlo = np.array([float((x.mid()
                                           - arb(float(x.mid()))).mid())
                                    for x in r['D']])
                    ldlt_arrays[f'{key}_D{sfx}_hi'] = dhi
                    ldlt_arrays[f'{key}_D{sfx}_lo'] = dlo
        if not ok:
            break
    rho_cert = bool(rho_runs.get('0.34', {}).get('certified'))
    tightest = None
    for t in ladder:
        if rho_runs.get(str(t), {}).get('certified'):
            tightest = t
    if ldlt_arrays:
        np.savez_compressed(LDLT_OUT, **ldlt_arrays)
    npz_check = check_ldlt_dump(At, R_at)
    print('ldlt dump check:', json.dumps(npz_check, default=str)[:300],
          flush=True)
    # Headline verdicts rest on the dump check (it folds the measured
    # factorisation residual S = M - L D L^T into the perturbation
    # bound); the ladder values are the search record.
    rho_cert = bool(npz_check.get('cat', {}).get('ok'))
    if npz_check.get('tight', {}).get('ok'):
        tightest = npz_check['tight']['t']
    elif rho_cert:
        tightest = 0.34
    else:
        tightest = None
    # (H-A) tail
    tail_lo, parts = tail_profile_lower(float(g.lower()), ball_up(qd),
                                        ball_up(asm['Cz']), asm['invF_up'],
                                        asm['uY_up'])
    cA_cert = min(cA_fin, tail_lo)
    res = {
        'date': '2026-11-15',
        'run': 'flint run 16: (H-A)/(H-rho) certification',
        'validations': val,
        'q_delta': str(qd),
        'g': str(g),
        'C_z': str(asm['Cz']),
        'invF_up': asm['invF_up'],
        'uY_up': asm['uY_up'],
        'HA_finite': {'cA_lower': cA_fin, 'argmin_j': jmin,
                      'profile_float_min_midpoints': prof_float_min,
                      'run13_float_value': 1.7372355768219303},
        'HA_tail': {'profile_lower_all_j_gt_812': tail_lo, 'parts': parts},
        'HA': {'catalogued_cA': 1.7372355768219303,
               'certified_cA': cA_cert,
               'certified_as_stated': bool(cA_cert >= 1.7372355768219303),
               'binding': 'finite' if cA_fin <= tail_lo else 'tail'},
        'Hrho': {'certified_rho_le_0.34': rho_cert,
                 'tightest_certified_threshold': tightest,
                 'ladder': rho_runs,
                 'ldlt_dump_check': npz_check,
                 'rho_float_midpoints': rho_float,
                 'run13_float_value': 0.33419022381006247,
                 'run13_value_status': 'SUPERSEDED -- float64 artifact '
                 '(see run note section 4): Arb-direct u_j balls exclude '
                 'run 13 u_all at j = 73, 74, 200, 412; the D-normalised '
                 'top form is a 1.14e6 cancellation in M - uu^T/g, and '
                 'the certified bound is the ladder above.'},
        'elapsed_s': time.time() - T00,
    }
    if not from_dumps:
        json.dump(res, open(RESULTS_OUT, 'w'), indent=1, default=str)
        dump = {'Q': {str(k): ser(v) for k, v in Q.items()},
                'A_diag': {str(j): ser(A[j - Z0, j - Z0])
                           for j in range(Z0, JMAX + 1)},
                'pivots_1.34I_minus_Atilde': [],
                'pivots_Atilde_minus_0.66I': [],
                'rho_ladder': rho_runs,
                'profile_lower': {str(k): v for k, v in prof_lo.items()},
                'g': ser(g), 'q_delta': ser(qd), 'C_z': ser(asm['Cz']),
                'invF_up': asm['invF_up'], 'uY_up': asm['uY_up']}
        json.dump(dump, open(BALLS_OUT, 'w'))
    print(json.dumps(res, indent=1, default=str))
    return res


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'certify'
    if mode == 'k2':
        compute_k2_all()
    elif mode == 'certify':
        run_certification(from_dumps=False)
    elif mode == 'verify':
        res = run_certification(from_dumps=True)
        ref = json.load(open(RESULTS_OUT))
        checks = {
            'HA_finite_cA_matches': abs(res['HA_finite']['cA_lower']
                                       - ref['HA_finite']['cA_lower']) < 1e-9,
            'HA_tail_matches': abs(res['HA_tail']['profile_lower_all_j_gt_812']
                                   - ref['HA_tail']['profile_lower_all_j_gt_812']) < 1e-9,
            'Hrho_matches': res['Hrho']['certified_rho_le_0.34']
            == ref['Hrho']['certified_rho_le_0.34'],
            'Hrho_tightest_matches': res['Hrho']
            ['tightest_certified_threshold']
            == ref['Hrho']['tightest_certified_threshold'],
            'ldlt_dump_ok': all(v['ok'] for v in
                                res['Hrho']['ldlt_dump_check'].values()),
            'validations_ok': res['validations']['V_all_ok'],
        }
        print('VERIFY:', json.dumps(checks))
        if not all(checks.values()):
            raise SystemExit('VERIFY FAILED')
    else:
        raise SystemExit(f'unknown mode {mode}')


if __name__ == '__main__':
    main()
