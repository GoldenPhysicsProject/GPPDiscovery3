# Flint, 2026-11-11: numerics-arm audit of Astra's completed-kernel factor-map note
# (authoritative copy: Astra_Completed_Kernel_Factor_Map_Obstruction_2026-11-11.md
#  in the RH goal's smatrix-lane folder).
# All tests use REAL test functions f. Conventions exactly as in the note:
#   q_{f,g}(u) = int f(v) conj(g(v-u)) dv ;  F(t) = int f(u) e^{-itu} du
#   M_f(s) = int f(u) e^{(s-1/2)u} du ; L+ = M_f(1), L- = M_f(0)
#   m_inf(t) = Re psi(1/4 + it/2) - log pi ; a_{p,k} = k log p ; w_{p,k} = (log p) p^{-k/2}
#   prime channel k: (P f)_k = sqrt(w/2)[f(u)-f(u+a)], (C f)_k = sqrt(w/2)[f(u)+f(u+a)]
#                    on I_{R,a} = [-R, R-a]
#   archimedean: P_inf f = (sqrt(m+) F, Lc), C_inf f = (sqrt(m-) F, Ls),
#                Lc = (L+ + L-)/sqrt2, Ls = (L+ - L-)/sqrt2
# Quadrature: uniform grid in u, trapezoid (f vanishes at grid ends);
#             uniform grid in t, trapezoid; psi via mpmath.
# Label: NUMERICAL (float64 grids; internal cross-checks printed).
import numpy as np
from mpmath import mp, digamma as mp_digamma

mp.dps = 30
LOGPI = np.log(np.pi)

def bump(u, c, w):
    x = (u - c) / w
    out = np.zeros_like(u)
    m = np.abs(x) < 1.0
    out[m] = np.exp(-1.0 / (1.0 - x[m] ** 2))
    return out

def m_inf(t):
    # Re psi(1/4 + i t/2) - log pi, evaluated pointwise with mpmath
    vals = np.empty_like(t)
    for i, ti in enumerate(t):
        vals[i] = float(mp_digamma(mp.mpc(0.25, ti / 2.0)).real) - LOGPI
    return vals

def fourier(f, u, t):
    du = u[1] - u[0]
    F = np.empty(len(t), dtype=complex)
    for j in range(0, len(t), 400):
        tt = t[j:j + 400]
        E = np.exp(-1j * np.outer(tt, u))
        F[j:j + 400] = (E @ f) * du
    return F

def primes_upto(n):
    sieve = np.ones(int(n) + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return np.nonzero(sieve)[0]

def shifted(f, u, a):
    return np.interp(u + a, u, f, left=0.0, right=0.0)

def analyse(f, u, R, t, mv, label):
    """Full carrier computation for one real test f supported in [-R,R]."""
    du = u[1] - u[0]
    dt = t[1] - t[0]
    Lp = np.trapezoid(f * np.exp(u / 2.0), u)
    Lm = np.trapezoid(f * np.exp(-u / 2.0), u)
    Lc = (Lp + Lm) / np.sqrt(2.0)
    Ls = (Lp - Lm) / np.sqrt(2.0)
    F = fourier(f, u, t)
    F2 = np.abs(F) ** 2
    mp_ = np.maximum(mv, 0.0)
    mm_ = np.maximum(-mv, 0.0)
    archP = np.trapezoid(mp_ * F2, t) / (2 * np.pi)
    archC = np.trapezoid(mm_ * F2, t) / (2 * np.pi)
    arch_signed = np.trapezoid(mv * F2, t) / (2 * np.pi)
    poleP, poleC = Lc ** 2, Ls ** 2
    pole_diff = poleP - poleC
    pole_Q = 2.0 * Lp * Lm  # L+ conj(L-) + L- conj(L+) for real f
    primeP = primeC = primeQ = 0.0
    channels = []
    for p in primes_upto(np.exp(2 * R)):
        lp = np.log(p)
        k = 1
        while k * lp <= 2 * R + 1e-12:
            a = k * lp
            w = lp * p ** (-k / 2.0)
            fp = shifted(f, u, a)    # f(u+a)
            fm = shifted(f, u, -a)   # f(u-a)
            q_a = np.trapezoid(f * fm, u)    # q(a)  = int f(v) f(v-a)
            q_ma = np.trapezoid(f * fp, u)   # q(-a) = int f(v) f(v+a)
            mask = (u >= -R - 1e-12) & (u <= R - a + 1e-9)
            us = u[mask]
            fs = f[mask]
            fps = fp[mask]
            cP = (w / 2.0) * np.trapezoid((fs - fps) ** 2, us)
            cC = (w / 2.0) * np.trapezoid((fs + fps) ** 2, us)
            primeP += cP
            primeC += cC
            primeQ += -w * (q_a + q_ma)
            channels.append((int(p), k, a, w, cP - cC, -w * (q_a + q_ma), q_a, q_ma))
            k += 1
    Pn2 = poleP + archP + primeP
    Cn2 = poleC + archC + primeC
    Q_direct = pole_Q + arch_signed + primeQ
    res = dict(label=label, R=R, Pn2=Pn2, Cn2=Cn2, diff=Pn2 - Cn2,
               Q=Q_direct, ratio=np.sqrt(Cn2 / Pn2),
               poleP=poleP, poleC=poleC, archP=archP, archC=archC,
               primeP=primeP, primeC=primeC, arch_signed=arch_signed,
               pole_Q=pole_Q, primeQ=primeQ, channels=channels)
    print(f"[{label}] R={R}: ||P||^2={Pn2:.10f} ||C||^2={Cn2:.10f} "
          f"diff={Pn2-Cn2:.10f} Q_direct={Q_direct:.10f} "
          f"|diff-Q|={abs(Pn2-Cn2-Q_direct):.3e} ratio={np.sqrt(Cn2/Pn2):.8f}")
    print(f"    pieces P: pole={poleP:.8f} arch={archP:.8f} prime={primeP:.8f} | "
          f"C: pole={poleC:.8f} arch={archC:.8f} prime={primeC:.8f}")
    return res

def main():
    out = []
    # ---------- (i) prime-channel identity, concrete bump ----------
    print("=== (i) prime-channel identity: f = bump(c=0.3, w=0.4), R = 2 ===")
    du = 0.002
    u = np.arange(-2.0, 2.0 + du, du)
    f = bump(u, 0.3, 0.4)
    t = np.arange(-100.0, 100.0 + 0.02, 0.02)
    mv = m_inf(t)
    r = analyse(f, u, 2.0, t, mv, "i-bump")
    for (p, k, a, w, lhs, rhs, qa, qma) in r["channels"]:
        if (p, k) in [(2, 1), (2, 2), (3, 1), (5, 1), (2, 3)]:
            alt_plus = +w * (qa + qma)
            alt_w2 = -w * w * (qa + qma)
            rel = f"{abs(lhs-rhs)/abs(rhs):.2e}" if rhs != 0 else "rhs=0"
            print(f"  p={p} k={k} a={a:.6f} w={w:.10f} LHS={lhs:.10f} "
                  f"-w(q+q-)={rhs:.10f} relerr={rel} "
                  f"+w(q+q-)={alt_plus:.10f} -w^2(q+q-)={alt_w2:.10f}")

    # ---------- (iii) archimedean pair ----------
    print("=== (iii) archimedean: f = bump(c=0, w=1), R = 2 ===")
    f1 = bump(u, 0.0, 1.0)
    r1 = analyse(f1, u, 2.0, t, mv, "iii-bump")
    for (p, k, a, w, lhs, rhs, qa, qma) in r1["channels"]:
        if (p, k) in [(2, 1), (2, 2), (3, 1), (5, 1), (7, 1)]:
            rel = f"{abs(lhs-rhs)/abs(rhs):.2e}" if rhs != 0 else "rhs=0"
            print(f"  wide-bump channel p={p} k={k} a={a:.6f} w={w:.10f} LHS={lhs:.10f} "
                  f"-w(q+q-)={rhs:.10f} relerr={rel}")
    gamma_diff = r1["archP"] - r1["archC"]
    print(f"  arch diff (P-C) = {gamma_diff:.10f}; signed int m|F|^2/2pi = "
          f"{r1['arch_signed']:.10f}; pole diff = {r1['poleP']-r1['poleC']:.10f} "
          f"vs pole Q = {r1['pole_Q']:.10f}")
    print(f"  collector-reading gamma diff (multipliers m and |m|, not sqrt): "
          f"int (|m|^2 - m^2)|F|^2/2pi = 0 exactly (m real); swapped-sqrt reading "
          f"gives {-r1['arch_signed']:.10f}")
    print(f"  m_inf(0) = {float(mp_digamma(mp.mpf('0.25')).real) - LOGPI:.10f}")

    # ---------- (ii) Theorem 5 window: R = 0.45, log2 < 2R < log3 ----------
    print("=== (ii) Theorem 5: R=0.45, a=log2, b = bump(c=-0.35, w=0.08) ===")
    R5 = 0.45
    a2 = np.log(2.0)
    w21 = a2 * 2.0 ** (-0.5)
    du5 = 0.0005
    u5 = np.arange(-R5, R5 + du5, du5)
    b = bump(u5, -0.35, 0.08)
    fb = b + shifted(b, u5, -a2)  # b(u) + b(u-a)
    mask = (u5 >= -R5 - 1e-12) & (u5 <= R5 - a2 + 1e-9)
    us = u5[mask]
    fbs = fb[mask]
    # f_b evaluated exactly at us + a (no interpolation): f_b(x) = b(x) + b(x - a),
    # b the bump centred -0.35, so f_b(x) = bump(x;-0.35) + bump(x;-0.35+a).
    fps = bump(us + a2, -0.35, 0.08) + bump(us + a2, -0.35 + a2, 0.08)
    Pcomp = np.sqrt(w21 / 2.0) * (fbs - fps)
    Ccomp = np.sqrt(w21 / 2.0) * (fbs + fps)
    pred = np.sqrt(2.0 * w21) * b[mask]
    print(f"  w_{{2,1}} = {w21:.12f}; sqrt(2 w) = {np.sqrt(2*w21):.12f}; "
          f"2^(1/4)*sqrt(log2) = {2**0.25*np.sqrt(a2):.12f}")
    print(f"  max|P_{{2,R}} f_b| on I_{{R,a}} = {np.max(np.abs(Pcomp)):.3e}")
    print(f"  max|C_{{2,R}} f_b - sqrt(2w) b| = {np.max(np.abs(Ccomp - pred)):.3e}")
    print(f"  ||b||^2 = {np.trapezoid(b**2, u5):.10f}; "
          f"||C_{{2,R}} f_b||^2 = {np.trapezoid(Ccomp**2, us):.10f}; "
          f"2 w ||b||^2 = {2*w21*np.trapezoid(b**2, u5):.10f}")
    print(f"  collector coefficient (2*w_coll, w_coll=(log2)/2) = {a2:.10f}; "
          f"ratio exact/collector = {np.sqrt(2*w21)/a2:.8f}")
    # full carrier for the same f_b (CGC in the Theorem-5 window)
    t5 = np.arange(-200.0, 200.0 + 0.05, 0.05)
    mv5 = m_inf(t5)
    analyse(fb, u5, R5, t5, mv5, "ii-fb-full")
    # moment-constrained b (Theorem 5 restriction): three bumps, null vector
    # of the 2x3 moment map (int b e^{u/2}, int b e^{-u/2}) via SVD
    B1 = bump(u5, -0.40, 0.045)
    B2 = bump(u5, -0.31, 0.045)
    B3 = bump(u5, -0.355, 0.03)
    M3 = np.array([[np.trapezoid(B1 * np.exp(u5 / 2), u5), np.trapezoid(B2 * np.exp(u5 / 2), u5),
                    np.trapezoid(B3 * np.exp(u5 / 2), u5)],
                   [np.trapezoid(B1 * np.exp(-u5 / 2), u5), np.trapezoid(B2 * np.exp(-u5 / 2), u5),
                    np.trapezoid(B3 * np.exp(-u5 / 2), u5)]])
    # null vector via SVD
    _, _, Vt = np.linalg.svd(M3)
    c3 = Vt[-1]
    bc = c3[0] * B1 + c3[1] * B2 + c3[2] * B3
    print(f"  3-bump constrained b: moments = "
          f"{np.trapezoid(bc*np.exp(u5/2), u5):.3e}, {np.trapezoid(bc*np.exp(-u5/2), u5):.3e}")
    fbc = bc + shifted(bc, u5, -a2)
    analyse(fbc, u5, R5, t5, mv5, "ii-fb-constrained")

    # ---------- (iv) CGC probe family ----------
    print("=== (iv) CGC probe ===")
    for (cc, ww, RR) in [(0.0, 1.0, 1.25), (0.0, 1.0, 2.0), (0.0, 1.0, 3.0), (0.0, 1.0, 4.0),
                         (-0.8, 0.5, 2.0), (-0.8, 0.5, 3.0), (1.2, 0.6, 2.0), (1.2, 0.6, 3.0),
                         (0.0, 0.35, 2.0), (0.0, 1.8, 3.0)]:
        uu = np.arange(-RR, RR + du, du)
        ff = bump(uu, cc, ww)
        analyse(ff, uu, RR, t, mv, f"iv-bump c={cc} w={ww}")
    # autocorr-type datum: f = q_{g,g} for g = bump(0, 0.7), support (-1.4, 1.4)
    for RR in (2.0, 3.0, 4.0):
        uu = np.arange(-RR, RR + du, du)
        g = bump(uu, 0.0, 0.7)
        # q(u) = int g(v) g(v-u) dv via FFT correlation on the grid
        n = len(g)
        G = np.fft.rfft(g, 2 * n)
        corr = np.fft.irfft(G * np.conj(G), 2 * n)[:n] * du  # corr[lagstep], lag>=0
        pos = uu >= 0
        lagsteps = np.round(uu[pos] / du).astype(int)
        fac = np.interp(np.abs(uu), uu[pos], corr[lagsteps])
        analyse(fac, uu, RR, t, mv, "iv-autocorr")

    # ---------- refinement check on the headline case ----------
    print("=== refinement: f=bump(0,1), R=2, du=0.001, dt=0.01, tmax=140 ===")
    du2 = 0.001
    u2 = np.arange(-2.0, 2.0 + du2, du2)
    t2 = np.arange(-140.0, 140.0 + 0.01, 0.01)
    mv2 = m_inf(t2)
    analyse(bump(u2, 0.0, 1.0), u2, 2.0, t2, mv2, "refine")

if __name__ == "__main__":
    main()
