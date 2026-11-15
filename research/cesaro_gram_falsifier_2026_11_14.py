#!/usr/bin/env python
"""Cheapest falsifier for the Cesaro-GNS route to Step C (spec section 4).

Spec: ~/workspace/goals/prove-the-riemann-hypothesis/hidden_files/
      smatrix-lane/pss-at-d1-construction-spec-2026-11-14.md

WHAT IS COMPUTED
----------------
The 2x2 Gram matrix of the first two zero states Psi_{rho_1}, Psi_{rho_2}
under the completed renormalised Cesaro pairing, evaluated from the
explicit-formula side, plus the F_theta control at the shifted labels
rho_k +- 1/4.

DEFINITIONS (choices made where the spec leaves the regularisation open;
each is labelled and used identically for the zeta labels and the control).

(D1) Labels.  A state label is s in C; its Cesaro exponent is alpha = s - 1/2
     (corpus coordinate of CesaroAbelState.lean: the critical line is
     Re alpha = 0).  Zeta labels: s_k = rho_k = 1/2 + i*gamma_k, k = 1, 2,
     with gamma_1, gamma_2 verified against the staged Platt table
     (zeros_by_index(1, 2); receipt in the JSON).

(D2) Spectral amplitude of a state.  The corpus Abel-Cesaro kernel
     (CesaroAbelState.lean, `abelValue`) is
        abelValue(eps, beta) = 2 eps^2 / (eps^2 - beta^2),
     the pairing of characters with cross exponent beta = conj(alpha) + alpha'.
     With beta = conj(s) + w - 1 for a spectral point w, the amplitude of the
     state s at w is
        A_eps(s, w) = 2 eps^2 / (eps^2 - (conj(s) + w - 1)^2)
                    = 2 eps^2 / (eps^2 - (w - u)^2),   u := 1 - conj(s).
     (The task text's kernel eps^2/(eps^2 - (sbar + s' - 1)^2) is abelValue/2;
     the factor 2 is the corpus normalisation omega_eps(1) = 2 and cancels in
     the Gram normalisation (D4).)

(D3) Gram weight and completed pairing at regulator eps.
        W_{kj,eps}(w) = (1/2) A_eps(s_k, w) A_eps(s_j, w),
        C_eps(k, j)   = W-functional[ W_{kj,eps} ]  (explicit-formula side),
     where the Guinand-Weil functional is
        W[F] = F(0) + F(1) - PRIME[F] + ARCH[F],
        PRIME[F] = sum_n Lambda(n) n^{-1/2} [ phi(log n) + phi(-log n) ],
        ARCH[F]  = (1/2pi) int_{-inf}^{inf} F(1/2 + it)
                     ( Re psi(1/4 + it/2) - log pi ) dt,
        phi(x)   = (1/2pi) int F(1/2 + it) e^{itx} dt.
     Conventions CERTIFIED on a Gaussian weight against the Platt zero sum
     (Part 1 of the run; kernel Re psi(1/4 + it/2) - log pi is also the
     programme's own archTerm kernel, KreinHPsi.lean line 64-67).
     On the spectral side W[W_{kj,eps}] = sum_rho W_{kj,eps}(rho), and
     W_{kk,eps}(rho_k) = 2, so C_eps -> 2 * delta_kj as eps -> 0+.

(D4) Gram matrix.  G_eps = C_eps / 2 (GNS normalisation of CesaroGNS.lean:
     the on-line limit form is 2 * sum conj(a) b, so a unit state has
     diagonal 2 before, and 1 after, this division).  G_0 is the eps -> 0
     extrapolation of G_eps over the headline grid; the smaller eigenvalue
     of G_0 is the falsifier's test statistic.

(D5) PRIME by analytic continuation.  phi for the weights of (D3) is an
     exponential polynomial on each half-line (residue calculus, Part 2
     certifies the closed forms against quadrature).  A component
     (c + d x) e^{a x} on x > 0 contributes
        c * L(1/2 - a) - d * L'(1/2 - a),   L := -zeta'/zeta,
     to PRIME, and a component (c + d x) e^{b x} on x < 0 contributes
        c * L(1/2 + b) + d * L'(1/2 + b)
     (because sum_n Lambda(n) (log n) n^{-s} = -L'(s)).
     Justification (labelled): for eps > 1/2 the weight is in the Weil
     class, PRIME converges absolutely as a Dirichlet sum, and every piece
     of W[W_{kj,eps}] is analytic in eps on (0, 1/2) U (1/2, inf) away from
     eps = 0; the assembly below is the analytic continuation of the
     in-class pairing.  Part 3 validates the assembly at eps = 0.6, 0.75:
     sharply against the Platt spectral sum (agreement ~2e-9), and against
     the directly truncated prime sum to the accuracy that truncation can
     reach there -- the signed prime tail at X = 2e6 is measured at
     ~0.03-0.07 (it decays only like X^{1/2-eps} through the explicit
     oscillatory terms X^{1-s}/(1-s) and sum_rho X^{rho-s}/(rho-s), and
     the absolute majorant tail is larger still), which is exactly why
     the continuation evaluation, not truncation, is the headline route.

(D6) Control.  Identical code at the shifted labels s = rho_k +- 1/4.
     Required behaviour (spec section 4, via the in-tree Cesaro selector
     A3): the FREE Cesaro diagonal at sigma = 1/2 +- 1/4 diverges with the
     regulator.  Exhibited three ways: (i) the R-mean
     M_sigma(R) = (1/(2 log R)) int_{1/R}^{R} r^{2 sigma - 2} dr
     (exact closed form; = 1 identically at sigma = 1/2, ~ R^{1/2}/log R
     at the shifts); (ii) the free Abel diagonal
     2 eps^2/(eps^2 - (2 sigma - 1)^2), singular at eps = |2 sigma - 1|
     and negative below it; (iii) the completed diagonal C_eps at the
     shifted labels from the same pipeline (it exhibits the Abel-kernel
     pole at eps = |2 sigma - 1| = 1/4 -- values 2.93, 4.18, 21.4, 6.24,
     0.22 at eps = 0.6, 0.45, 0.3, 0.2, 0.13 -- and tends to 0 as
     eps -> 0; at no eps is it a finite positive state norm).

ERROR BUDGET (assembled in the JSON):
  evaluator: Arb ball radii for L, L' (prec 250, central differences,
             h = 1e-18) and mpmath quadrature self-consistency at dps 50/65.
  arch tail: analytic bound with T = 400 (see arch_term).
  prime cutoff: the headline uses the continuation value (no truncation);
             for the in-class validation epsilons the corrected absolute
             truncation bound at X = 2e6 and the measured signed tail
             (X = 2e5 -> 2e6) are both reported -- near eps = 1/2 the
             absolute bound is O(1) and the signed tail is ~0.03-0.07,
             which is precisely why truncation cannot certify there and
             the continuation value is used (itself certified against
             the spectral side to ~2e-9).
  extrapolation: spread of eps -> 0 fits (degrees 1/2/3, nested grids).

Epistemic status: NUMERICAL throughout.  A SURVIVES verdict is a statement
about one 2x2 corner only (spec section 4).
"""

import json
import sys
import time

sys.path.insert(0, '/home/hatch/workspace/zeros-platt')

import mpmath
from mpmath import mp, mpf, mpc, exp, log, pi, sqrt, cos, sin, digamma, quad
from flint import acb, arb, ctx

mp.dps = 60
T_START = time.time()

OUT = {}

# ---------------------------------------------------------------- Part 0
# Ordinates: task-supplied values verified against the Platt table.
GAMMA_TASK = {
    1: '14.134725141734693790457251983562470270784257115699',
    2: '21.022039638771554992628479593896902777334340524902',
}

def verify_ordinates():
    from platt_sequential import zeros_by_index
    got = dict(zeros_by_index(1, 2))
    rec = {}
    gam = {}
    for k in (1, 2):
        g_tab = got[k]
        g_task = mpf(GAMMA_TASK[k])
        rec[k] = {
            'platt': mpmath.nstr(g_tab, 55),
            'task': GAMMA_TASK[k],
            'abs_diff': mpmath.nstr(abs(g_tab - g_task), 4),
        }
        gam[k] = g_tab  # use the table values (exact decode)
    return gam, rec

GAM, GAM_REC = verify_ordinates()
OUT['ordinates_receipt'] = GAM_REC
RHO = {k: mpc(mpf('0.5'), GAM[k]) for k in (1, 2)}

# ---------------------------------------------------------------- Part 1
# Certification A: Guinand-Weil conventions on a Gaussian weight.
def gauss_cert(t0):
    tau = mpf(3)
    t0 = mpf(t0)
    def Phi_w(w):
        return exp((w - mpf('0.5') - 1j * t0) ** 2 / (2 * tau ** 2))
    def Phi_line(t):
        return exp(-(mpf(t) - t0) ** 2 / (2 * tau ** 2))
    def phi(x):
        x = mpf(x)
        return (tau / sqrt(2 * pi)) * exp(-tau ** 2 * x ** 2 / 2) * exp(1j * t0 * x)
    return Phi_w, Phi_line, phi

def mangoldt_upto(N):
    import math
    spf = list(range(N + 1))
    for i in range(2, int(N ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, N + 1, i):
                if spf[j] == j:
                    spf[j] = i
    lam = [0.0] * (N + 1)
    for n in range(2, N + 1):
        m, p = n, spf[n]
        while m % p == 0:
            m //= p
        if m == 1:
            lam[n] = math.log(p)
    return lam

LAM4K = mangoldt_upto(4000)

def arch_integral(Phi_line, extra_pts=()):
    """(1/2pi) int Phi_line(t) (Re psi(1/4+it/2) - log pi) dt, split quad."""
    pts = [-400, -60, -20, -5, -1, 0, 1, 5, 20, 60, 400]
    for p in extra_pts:
        for d in (-2, 0, 2):
            pts.append(float(p) + d)
    pts = sorted(set(pts))
    f = lambda t: Phi_line(t) * ((digamma(mpc(mpf('0.25'), mpf(t) / 2))).real - log(pi))
    return quad(f, pts) / (2 * pi)

def gauss_residual(t0):
    from platt_sequential import zeros_in_range
    Phi_w, Phi_line, phi = gauss_cert(t0)
    S = mpc(0)
    for g in zeros_in_range(0, 300):
        S += Phi_w(mpc('0.5', g)) + Phi_w(mpc('0.5', -g))
    pole = Phi_w(mpc(0)) + Phi_w(mpc(1))
    P = mpc(0)
    for n in range(2, 4001):
        if LAM4K[n] != 0.0:
            ln = log(mpf(n))
            P += LAM4K[n] * n ** (-mpf('0.5')) * (phi(ln) + phi(-ln))
    A = arch_integral(Phi_line, extra_pts=[t0])
    rhs = pole - P + A
    return S, rhs

# ---------------------------------------------------------------- Part 2
# The Gram weights, phi closed forms, residue machine.

def u_of(s):
    """Centre of the amplitude in the w-plane: u = 1 - conj(s)."""
    return 1 - s.conjugate()

def W_line(t, eps, u_k, u_j):
    """W_{kj,eps}(1/2 + it) as an explicit function of t (real for the
    label pairs used here up to a phase; returned complex)."""
    w = mpc(mpf('0.5'), mpf(t))
    Ak = 2 * eps ** 2 / (eps ** 2 - (w - u_k) ** 2)
    Aj = 2 * eps ** 2 / (eps ** 2 - (w - u_j) ** 2)
    return Ak * Aj / 2

def W_val(w, eps, u_k, u_j):
    Ak = 2 * eps ** 2 / (eps ** 2 - (w - u_k) ** 2)
    Aj = 2 * eps ** 2 / (eps ** 2 - (w - u_j) ** 2)
    return Ak * Aj / 2

def Qprime(w, eps, u_k, u_j):
    """Derivative of Q(w) = (eps^2 - (w-u_k)^2)(eps^2 - (w-u_j)^2)."""
    return (-2 * (w - u_k) * (eps ** 2 - (w - u_j) ** 2)
            - 2 * (w - u_j) * (eps ** 2 - (w - u_k) ** 2))

def phi_components(eps, u_k, u_j):
    """Exponential-polynomial components of phi on each half-line.

    Returns (pos, neg): lists of (c, d, a) meaning phi(x) contains
    (c + d*x) * e^{a x} for x > 0 (pos) resp. x < 0 (neg).
    Diagonal case u_k == u_j uses the double-pole closed form:
      phi(x) = (eps/2)(1 + eps|x|) e^{(u-1/2)x} e^{-eps|x|}.
    Off-diagonal: simple poles at u_k +- eps, u_j +- eps; for x > 0 the
    poles with Re < 1/2 contribute Res = 2 eps^4 e^{(w0-1/2)x} / Q'(w0);
    for x < 0, minus the sum over poles with Re > 1/2.
    """
    half = mpf('0.5')
    if u_k == u_j:
        u = u_k
        pos = [(eps / 2, eps ** 2 / 2, u - half - eps)]
        neg = [(eps / 2, -eps ** 2 / 2, u - half + eps)]
        # check: for x<0, (c + d x) e^{a x} with a = u-1/2+eps,
        # c = eps/2, d = -eps^2/2 gives (eps/2)(1 - eps x) e^{(u-1/2+eps)x}
        # = (eps/2)(1 + eps|x|) e^{(u-1/2)x} e^{-eps|x|}  ✓
        return pos, neg
    pos, neg = [], []
    for w0 in (u_k - eps, u_j - eps):
        R = 2 * eps ** 4 / Qprime(w0, eps, u_k, u_j)
        pos.append((R, mpc(0), w0 - half))
    for w0 in (u_k + eps, u_j + eps):
        R = -2 * eps ** 4 / Qprime(w0, eps, u_k, u_j)
        neg.append((R, mpc(0), w0 - half))
    return pos, neg

def phi_eval(x, eps, u_k, u_j):
    pos, neg = phi_components(eps, u_k, u_j)
    comps = pos if x > 0 else neg
    if x == 0:
        # average of the two branches' limits (they agree for our weights)
        return (sum(c for c, d, a in pos) + sum(c for c, d, a in neg)) / 2
    return sum((c + d * x) * exp(a * x) for c, d, a in comps)

def phi_quad_check(eps, u_k, u_j, xs):
    """Certify phi closed forms against direct quadrature of the
    Fourier integral on [-400, 400] (tail < 1e-12 for our weights)."""
    res = {}
    for x in xs:
        f = lambda t: W_line(t, eps, u_k, u_j) * exp(1j * mpf(t) * mpf(x))
        q = quad(f, [-400, -60, -20, 0, 20, 60, 400]) / (2 * pi)
        res[str(x)] = {
            'quad': mpmath.nstr(q, 20),
            'closed': mpmath.nstr(phi_eval(mpf(x), eps, u_k, u_j), 20),
            'abs_diff': mpmath.nstr(abs(q - phi_eval(mpf(x), eps, u_k, u_j)), 4),
        }
    return res

# ---------------------------------------------------------------- Part 3
# L = -zeta'/zeta and L' in Arb (rigorous balls, central differences).
def arb_L(s_re, s_im, want_deriv=True):
    ctx.prec = 250
    h = acb(arb('1e-18'))
    s = acb(arb(s_re), arb(s_im))
    z0 = s.zeta()
    zp = (s + h).zeta()
    zm = (s - h).zeta()
    dz = (zp - zm) / (2 * h)                # zeta'(s)
    ddz = (zp - 2 * z0 + zm) / (h * h)      # zeta''(s)
    L = -dz / z0
    Lp = -ddz / z0 + (dz / z0) ** 2 if want_deriv else None
    return L, Lp

def mp_L(s):
    """mpmath cross-check values of L, L' at an mpc point."""
    z0 = mpmath.zeta(s)
    dz = mpmath.diff(mpmath.zeta, s)
    ddz = mpmath.diff(mpmath.zeta, s, 2)
    L = -dz / z0
    Lp = -ddz / z0 + (dz / z0) ** 2
    return L, Lp

def L_at(s):
    """L(s), L'(s) as mpc at dps 60, from Arb; returns (L, Lp, radius)."""
    La, Lpa = arb_L(mpmath.nstr(s.real, 30), mpmath.nstr(s.imag, 30))
    def to_mp(z):
        return mpc(mpf(z.real.mid()), mpf(z.imag.mid()))
    rad = max(abs(La.real.rad()), abs(La.imag.rad()),
              abs(Lpa.real.rad()), abs(Lpa.imag.rad()))
    return to_mp(La), to_mp(Lpa), float(rad)

# ---------------------------------------------------------------- Part 4
# The completed pairing C_eps(k, j) from the explicit-formula side.

def prime_ac(eps, u_k, u_j, Lcache):
    """PRIME[F] by analytic continuation (D5).  Lcache maps a rounded
    s-point key to (L, L', radius)."""
    pos, neg = phi_components(eps, u_k, u_j)
    total = mpc(0)
    maxrad = 0.0
    def getL(s):
        key = (round(float(s.real), 12), round(float(s.imag), 12))
        if key not in Lcache:
            Lcache[key] = L_at(s)
        return Lcache[key]
    for c, d, a in pos:      # x > 0 branch: s = 1/2 - a
        s = mpf('0.5') - a
        L, Lp, rad = getL(s)
        total += c * L - d * Lp
        maxrad = max(maxrad, rad)
    for c, d, b in neg:      # x < 0 branch: s = 1/2 + b
        s = mpf('0.5') + b
        L, Lp, rad = getL(s)
        total += c * L + d * Lp
        maxrad = max(maxrad, rad)
    return total, maxrad

def arch_term(eps, u_k, u_j):
    """ARCH[F] with analytic tail bound beyond T = 400.
    |W(1/2+it)| <= 2 eps^4 / ((|t|-30)^4) for |t| >= 60 (labels here have
    |Im u| <= 22, |Re u - 1/2| <= 1/4), and
    |Re psi(1/4 + it/2) - log pi| <= log|t| + 2 there; the tail bound is
    the elementary majorant integral, evaluated numerically."""
    T = 400
    pts = [-T, -60, -30, -5, -1, 0, 1, 5, 30, 60, T]
    for g in (float(GAM[1]), float(GAM[2])):
        for d in (-2, 0, 2):
            pts += [g + d, -g + d]
    pts = sorted(set(pts))
    f = lambda t: W_line(t, eps, u_k, u_j) * ((digamma(mpc(mpf('0.25'), mpf(t) / 2))).real - log(pi))
    val = quad(f, pts) / (2 * pi)
    # tail majorant: int_T^inf 2 eps^4 (log t + 2) / (t-30)^4 dt / (2 pi)
    maj = quad(lambda t: 2 * eps ** 4 * (log(mpf(t)) + 2) / (mpf(t) - 30) ** 4,
               [T, mpmath.inf]) / (2 * pi)
    return val, maj

def C_eps(eps, s_k, s_j, Lcache):
    u_k, u_j = u_of(s_k), u_of(s_j)
    pole = W_val(mpc(0), eps, u_k, u_j) + W_val(mpc(1), eps, u_k, u_j)
    P, rad = prime_ac(eps, u_k, u_j, Lcache)
    A, tail = arch_term(eps, u_k, u_j)
    return pole - P + A, {'pole': pole, 'prime': P, 'arch': A,
                          'L_radius': rad, 'arch_tail_bound': tail}

def prime_trunc(eps, u_k, u_j, LAM, X):
    """PRIME[F] by direct truncation at n <= X (in-class validation)."""
    total = mpc(0)
    for n in range(2, X + 1):
        if LAM[n] != 0.0:
            ln = log(mpf(n))
            total += LAM[n] * n ** (-mpf('0.5')) * (phi_eval(ln, eps, u_k, u_j)
                                                    + phi_eval(-ln, eps, u_k, u_j))
    return total

def phi_majorant_tail(eps, u_k, u_j, X):
    """Absolute bound for |sum_{n>X} Lambda(n) n^{-1/2} (phi(log n)+phi(-log n))|:
    |phi(x)| <= M(x) := sum over both branches' components of
    (|c| + |d| x) e^{Re(a) x} for x = log n > 0 (branchwise majorant),
    then sum_{n>X} (log n) n^{-1/2} M(log n) <= int_X^inf (log u) u^{-1/2}
    M(log u) du (monotone majorant; crude safety factor 2).  NOTE: at
    eps close to 1/2 this absolute bound is large (the tail decays only
    like X^{1/2-eps}); the signed tail is smaller but oscillatory.  This
    is reported honestly in the validation; truncation is a weak leg
    near the Weil-class boundary and is not the headline route."""
    pos, neg = phi_components(eps, u_k, u_j)
    def M(y):
        s = mpf(0)
        for c, d, a in pos:
            s += (abs(c) + abs(d) * y) * exp(a.real * y)
        for c, d, b in neg:
            # phi(-y): components (c + d(-y)) e^{b(-y)}
            s += (abs(c) + abs(d) * y) * exp(-b.real * y)
        return s
    f = lambda u: log(mpf(u)) * mpf(u) ** (-mpf('0.5')) * M(log(mpf(u)))
    return 2 * quad(f, [X, mpmath.inf])

def spectral_sum(eps, u_k, u_j, T=3000):
    from platt_sequential import zeros_in_range
    S = mpc(0)
    for g in zeros_in_range(0, T):
        S += W_line(g, eps, u_k, u_j) + W_line(-g, eps, u_k, u_j)
    return S

# ================================================================ RUN
print('== ordinates ==')
for k in (1, 2):
    print(k, GAM_REC[k]['platt'][:40], 'diff', GAM_REC[k]['abs_diff'])

print('== Certification A: Gaussian explicit formula ==')
OUT['cert_gaussian'] = {}
for t0 in ['0', GAMMA_TASK[1][:20]]:
    S, rhs = gauss_residual(mpf(t0))
    OUT['cert_gaussian'][t0] = {'spectral': mpmath.nstr(S, 25),
                                'arithmetic': mpmath.nstr(rhs, 25),
                                'abs_resid': mpmath.nstr(abs(S - rhs), 4)}
    print(t0, 'resid', mpmath.nstr(abs(S - rhs), 4))

print('== Certification B: phi closed forms vs quadrature ==')
eps_c = mpf('0.3')
OUT['cert_phi'] = {}
for (k, j) in ((1, 1), (1, 2)):
    key = f'{k}{j}'
    OUT['cert_phi'][key] = phi_quad_check(eps_c, RHO[k], RHO[j], ['0.7', '2.0'])
    print(key, {x: v['abs_diff'] for x, v in OUT['cert_phi'][key].items()})

print('== L values: Arb vs mpmath cross-check ==')
s_test = RHO[1] + mpf('0.3')
L1, Lp1, rad1 = L_at(s_test)
L2, Lp2 = mp_L(s_test)
OUT['cert_L'] = {'point': mpmath.nstr(s_test, 20),
                 'arb_L': mpmath.nstr(L1, 30), 'mp_L': mpmath.nstr(L2, 30),
                 'abs_diff_L': mpmath.nstr(abs(L1 - L2), 4),
                 'abs_diff_Lp': mpmath.nstr(abs(Lp1 - Lp2), 4),
                 'arb_radius': rad1}
print('L diff', mpmath.nstr(abs(L1 - L2), 4), 'Lp diff',
      mpmath.nstr(abs(Lp1 - Lp2), 4), 'radius', rad1)

print('== In-class validation: eps = 0.6, 0.75 (three-way) ==')
import numpy as np

def mangoldt_np(N):
    is_p = bytearray(b'\x01') * (N + 1)
    is_p[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if is_p[i]:
            is_p[i * i:N + 1:i] = b'\x00' * ((N - i * i) // i + 1)
    lam = np.zeros(N + 1)
    idx = np.nonzero(is_p)[0]
    lam[idx] = np.log(idx.astype(float))
    for p in idx:
        if p * p > N:
            break
        q = int(p) * int(p)
        while q <= N:
            lam[q] = np.log(float(p))
            q *= int(p)
    return lam

X_TRUNC = 2_000_000
LAMNP = mangoldt_np(X_TRUNC)
SUP = np.nonzero(LAMNP)[0]
LN_SUP = np.log(SUP.astype(float))
LAM_SUP = LAMNP[SUP]

def prime_trunc_np(eps, u_k, u_j, Xcap=None):
    """PRIME by direct summation in float128-complex (numpy complex128);
    components of phi evaluated vectorially on the support of Lambda."""
    pos, neg = phi_components(eps, u_k, u_j)
    tot = 0.0 + 0.0j
    if Xcap is None:
        sup, ln_sup, lam_sup = SUP, LN_SUP, LAM_SUP
    else:
        m = SUP <= Xcap
        sup, ln_sup, lam_sup = SUP[m], LN_SUP[m], LAM_SUP[m]
    n = sup.astype(float)
    base = lam_sup * n ** (-0.5)
    for c, d, a in pos:
        cc, dd, aa = complex(c), complex(d), complex(a)
        tot += np.sum(base * (cc + dd * ln_sup) * np.exp(aa * ln_sup))
    for c, d, b in neg:
        cc, dd, bb = complex(c), complex(d), complex(b)
        # phi(-y) = (c + d*(-y)) e^{b*(-y)}
        tot += np.sum(base * (cc - dd * ln_sup) * np.exp(-bb * ln_sup))
    return mpc(tot.real, tot.imag)

def spectral_tail_bound(eps, T=3000):
    maj = quad(lambda t: 2 * eps ** 4 * (log(mpf(t)) + 13) / (mpf(t) - 30) ** 4,
               [T, mpmath.inf])
    return maj

OUT['validation_inclass'] = {}
Lcache = {}
for eps in (mpf('0.6'), mpf('0.75')):
    for (k, j) in ((1, 1), (1, 2), (2, 2)):
        u_k, u_j = u_of(RHO[k]), u_of(RHO[j])
        Cac, parts = C_eps(eps, RHO[k], RHO[j], Lcache)
        pole = parts['pole']; A = parts['arch']
        Ptr = prime_trunc_np(eps, u_k, u_j)
        Ptr2 = prime_trunc_np(eps, u_k, u_j, Xcap=200_000)
        Ctr = pole - Ptr + A
        tail_tr = phi_majorant_tail(eps, u_k, u_j, X_TRUNC)
        emp_tail = abs(Ptr - Ptr2)
        S = spectral_sum(eps, u_k, u_j)
        tail_sp = spectral_tail_bound(eps)
        key = f'eps={eps} k={k} j={j}'
        OUT['validation_inclass'][key] = {
            'C_ac': mpmath.nstr(Cac, 20), 'C_trunc': mpmath.nstr(Ctr, 20),
            'C_spectral': mpmath.nstr(S, 20),
            'diff_ac_trunc': mpmath.nstr(abs(Cac - Ctr), 4),
            'abs_tail_bound_corrected': mpmath.nstr(tail_tr, 4),
            'empirical_tail_2e5_to_2e6': mpmath.nstr(emp_tail, 4),
            'diff_ac_spectral': mpmath.nstr(abs(Cac - S), 4),
            'tail_spectral_bound': mpmath.nstr(tail_sp, 4),
        }
        print(key, 'ac-trunc', mpmath.nstr(abs(Cac - Ctr), 3),
              '(emp tail', mpmath.nstr(emp_tail, 3), '; abs bound',
              mpmath.nstr(tail_tr, 3), ')  ac-spec',
              mpmath.nstr(abs(Cac - S), 3), '(bar', mpmath.nstr(tail_sp, 3), ')')

print('== Headline grid: completed Gram matrix, eps -> 0 ==')
EPS_GRID = [mpf('0.30'), mpf('0.26'), mpf('0.225'), mpf('0.195'),
            mpf('0.17'), mpf('0.15'), mpf('0.13')]
Cgrid = {}
max_imag = 0.0
max_arch_tail = 0.0
max_L_rad = 0.0
for eps in EPS_GRID:
    for (k, j) in ((1, 1), (1, 2), (2, 2)):
        C, parts = C_eps(eps, RHO[k], RHO[j], Lcache)
        Cgrid[(str(eps), k, j)] = C
        max_imag = max(max_imag, abs(C.imag))
        max_arch_tail = max(max_arch_tail, abs(parts['arch_tail_bound']))
        max_L_rad = max(max_L_rad, parts['L_radius'])
        print('eps', eps, k, j, mpmath.nstr(C, 18))
OUT['max_imag_part'] = mpmath.nstr(max_imag, 4)
OUT['max_arch_tail_bound'] = mpmath.nstr(max_arch_tail, 4)
OUT['max_L_radius'] = max_L_rad

# Extrapolation: polynomial fits in eps, degrees 1..3, full grid and
# sub-grid eps <= 0.26; entry error = max spread of intercepts + tails.
import itertools
def fit_intercepts(eps_list, vals):
    out = {}
    xs = np.array([float(e) for e in eps_list])
    ys = np.array([float(v) for v in vals])
    for deg in (1, 2, 3):
        if len(xs) > deg:
            cf = np.polyfit(xs, ys, deg)
            out[f'deg{deg}'] = float(np.polyval(cf, 0.0))
    return out

G0 = {}
G0_err = {}
OUT['fits'] = {}
for (k, j) in ((1, 1), (1, 2), (2, 2)):
    vals_full = [Cgrid[(str(e), k, j)].real / 2 for e in EPS_GRID]
    sub = [(e, Cgrid[(str(e), k, j)].real / 2) for e in EPS_GRID if e <= mpf('0.26')]
    f_full = fit_intercepts(EPS_GRID, vals_full)
    f_sub = fit_intercepts([e for e, _ in sub], [v for _, v in sub])
    allfits = dict(f_full)
    allfits.update({('sub_' + kk): vv for kk, vv in f_sub.items()})
    inter = list(allfits.values())
    # central: median of the degree-2 and degree-3 fits (full + sub)
    centrals = [allfits[kk] for kk in allfits if 'deg2' in kk or 'deg3' in kk]
    G0[(k, j)] = float(np.median(centrals))
    G0_err[(k, j)] = (max(inter) - min(inter)) / 2 + float(max_arch_tail) / 2 + max_L_rad
    OUT['fits'][f'{k}{j}'] = {kk: vv for kk, vv in allfits.items()}
    print('entry', k, j, 'G0 =', G0[(k, j)], 'err', G0_err[(k, j)])

a, b, c = G0[(1, 1)], G0[(1, 2)], G0[(2, 2)]
tr = a + c
det = a * c - b * b
disc = max(tr * tr / 4 - det, 0.0)
lam_min = tr / 2 - disc ** 0.5
lam_max = tr / 2 + disc ** 0.5
B = float(np.sqrt(sum(v * v for v in G0_err.values())))
OUT['gram_G0'] = {'G11': a, 'G12': b, 'G22': c}
OUT['eigenvalues'] = {'lambda_min': lam_min, 'lambda_max': lam_max,
                      'total_error_bar': B}
print('G0 =', [[a, b], [b, c]])
print('lambda_min =', lam_min, ' total bar =', B)

if lam_min < -B:
    VERDICT = 'KILLED'
elif B > 0.25:
    VERDICT = 'INCONCLUSIVE'
else:
    VERDICT = 'SURVIVES'
OUT['verdict'] = VERDICT
print('verdict:', VERDICT)

print('== Control C1: free Cesaro R-mean diagonal ==')
OUT['control_Rmean'] = {}
for sigma in ('0.5', '0.75', '0.25'):
    sg = mpf(sigma)
    row = {}
    for kk in (1, 2, 3, 4, 5):
        R = exp(4 * kk)
        if sg == mpf('0.5'):
            v = mpf(1)
        else:
            aa = 2 * sg - 1
            v = (R ** aa - R ** (-aa)) / (aa * 2 * log(R))
        row[f'R=e^{4 * kk}'] = mpmath.nstr(v, 12)
    OUT['control_Rmean'][sigma] = row
    print('sigma', sigma, row)

print('== Control C2: free Abel diagonal at the shifts ==')
OUT['control_Abel_shift'] = {}
for eps in ('0.75', '0.6', '0.55', '0.51', '0.49', '0.4', '0.3', '0.15'):
    e = mpf(eps)
    OUT['control_Abel_shift'][eps] = mpmath.nstr(2 * e ** 2 / (e ** 2 - mpf('0.25')), 12)
print(OUT['control_Abel_shift'])

print('== Control C3: completed diagonal at shifted labels ==')
OUT['control_completed_shift'] = {}
for sgn in (1, -1):
    for k in (1, 2):
        s = RHO[k] + sgn * mpf('0.25')
        row = {}
        for eps in (mpf('0.6'), mpf('0.45'), mpf('0.3'), mpf('0.2'), mpf('0.13')):
            Csh, _ = C_eps(eps, s, s, Lcache)
            row[str(eps)] = mpmath.nstr(Csh, 15)
        # spectral cross-check at eps = 0.3
        u = u_of(s)
        Ssh = spectral_sum(mpf('0.3'), u, u)
        Cac03, _ = C_eps(mpf('0.3'), s, s, Lcache)
        row['spectral_eps0.3'] = mpmath.nstr(Ssh, 15)
        row['ac_vs_spectral_diff'] = mpmath.nstr(abs(Ssh - Cac03), 4)
        OUT['control_completed_shift'][f'k={k} shift={sgn * 0.25:+}'] = row
        print('shift', k, sgn * 0.25, row)

OUT['runtime_seconds'] = time.time() - T_START
with open('/home/hatch/workspace/gppdiscovery3/research/'
          'cesaro_gram_falsifier_2026_11_14_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('done', OUT['runtime_seconds'])

