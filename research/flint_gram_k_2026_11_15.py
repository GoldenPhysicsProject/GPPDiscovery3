# Flint run 15, item 3' (2026-11-15): ball Gram(k) and CC^T via R10 identity CC^T = Gram(k) - HH^T,
# i.e. Gram(k) = HH^T + sum_{j>=73} C~_j C~_j^T,  C~_{a,j} = w_j q_a + om_j^2 rho_j(x_a)  (R8 (T2) split).
#  - H (j=1..72) and A_X' from Heart R7 Arb shards (GPPDiscovery5 1da974d, round7_turbo).
#  - q_{phi_i} = C~_{delta,i} = H[12][i-1] (identity q_x = <k_delta, phi_i>, verified vs spectral route);
#    q_delta from the kernel route q_delta = -1/2[K''(0)-K''(1/2)] (matches R8 float 0.1197039).
#  - rho via R12's Q-reduction: Q(w) = R(w) + Pi(w); R(om_j) by the exact finite-interval form
#    R(om_j) = (pi/2) int_0^{5/2} h(u) s(u) sin(om_j u) du (m-series + digamma; validated vs t2rig/R8).
#  - Head j=73..J0 (J0=4000) summed exactly in balls (blocked arb_mat products);
#    tail j>J0: rank-one part W*qq^T exact (W = sum_{j>J0} w_j^2), remainder bounded entrywise
#    via R12 envelopes: |om^4 rho_j(delta)| <= C_delta = 19.5771; |om^4 rho_j(phi_i)| <= 69.5492/(om_j-om_i)
#    (constants rounded UP from R12's Arb values C_delta = 19.57696116..., Cphi_num = 69.54894066...).
from flint import arb, acb, ctx, arb_mat
import json, time, subprocess, tarfile, io, os, glob
ctx.prec = 256
PI = arb.pi(); ELL = arb(1)/4; SQ8 = arb(8).sqrt()
BI2 = 16*(ELL).sinh()**2
SINH8 = (ELL/2).sinh()
NS = [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]
J0 = 4000
C_DELTA = arb('19.5771')
C_PHI = arb('69.5492')
SIGC = SQ8*SINH8  # |sig(j)| <= SIGC/om_j
def om(j): return 4*PI*j
def sig(j):
    w = om(j); return -SQ8*(-1)**j*w*SINH8/(arb(1)/4 + w*w)
def w_(j): return SQ8*(-1)**(j+1)/om(j)
# ---------------- R / Pi / Q ----------------
def E_(s, x, w): return (-s*x).exp()*(s*(w*x).sin() + w*(w*x).cos())/(s*s + w*w)
def I_piece(a, b, w, M=110):
    tot = arb(0)
    for m in range(M+1):
        s = arb(2*m) + arb(1)/2
        tot += 2*(E_(s, a, w) - E_(s, b, w))
    s0 = 2*(M+1) + arb(1)/2
    return tot + arb(0, 4*(-s0*a).exp()/s0/(1-(-2*a).exp()))
def I_first(w, M=110):
    z = acb(arb(1)/4) + acb(0,1)*(w/2)
    tot = z.digamma().imag
    for m in range(M+1):
        s = arb(2*m) + arb(1)/2
        tot -= 2*E_(s, arb(1)/2, w)
    s0 = 2*(M+1) + arb(1)/2
    return tot + arb(0, 4*(-s0/2).exp()/s0/(1-arb(-1).exp()))
PIECES = [(arb(0), arb(1)/2, arb(3)/16), (arb(1)/2, arb(1), arb(1)/8),
          (arb(1), arb(3)/2, -arb(1)/8), (arb(3)/2, arb(2), -arb(1)/32),
          (arb(2), arb(5)/2, arb(1)/32)]
def R_of(w):
    tot = arb(0)
    for (a,b,s) in PIECES:
        tot += s*(I_first(w) if a == 0 else I_piece(a, b, w))
    return (PI/2)*tot
CS = [arb(3)/16, -arb(1)/16, -arb(1)/4, arb(3)/32, arb(1)/16, -arb(1)/32]
def s_of(u):
    uf = float(u.mid()); tot = arb(0)
    for k in range(6):
        if k/2 < uf: tot += CS[k]
    return tot
def Pi_of(w):
    return PI*sum(((arb(p).log()/arb(n).sqrt())*s_of(arb(n).log())*(w*arb(n).log()).sin() for n, p in NS), arb(0))
def Q_of(w): return R_of(w) + Pi_of(w)
PREF = (16/PI)*SQ8
def w4rho_delta(j, Qj):
    return PREF*(-1)**j*Qj - (BI2**2/8)*SINH8*sig(j)
def w4rho_phi(i, j, Qi, Qj):
    wi, wj = om(i), om(j)
    return (128/PI)*(-1)**(i+j)*(wi*wj/(wj*wj - wi*wi))*(Qj/wj - Qi/wi) - (BI2**2/8)*sig(i)*sig(j)
# ---------------- K'' for q_delta (kernel route, eps-split integrals) ----------------
LOGN = {n: arb(n).log() for n, p in NS}
def beta2_val(y):
    a = abs(y); alo = float(a.lower()); ahi = float(a.upper())
    def f1(x): return -2 + 3*x
    def f2(x): return 2 - x
    if ahi < 1: return f1(a)
    if alo > 2: return arb(0)
    if alo > 1 and ahi < 2: return f2(a)
    cands = []
    if alo < 1: cands.append(f1(a))
    if ahi > 1 and alo < 2: cands.append(f2(a))
    if ahi > 2: cands.append(arb(0))
    lo = min(float(c.lower()) for c in cands); hi = max(float(c.upper()) for c in cands)
    return arb(lo) + arb(0, hi-lo) if hi > lo else cands[0]
def beta2_lin(c, sgn, umid):
    ymid = c + sgn*umid
    a = abs(float(ymid.mid()))
    if a > 2: return (arb(0), arb(0))
    if a < 1:
        if float((c + sgn*umid).mid()) >= 0: return (-2 + 3*c, 3*sgn)
        return (-2 - 3*c, -3*sgn)
    if float((c + sgn*umid).mid()) >= 0: return (2 - c, -sgn)
    return (2 + c, sgn)
MK = 130
def H0s(a, b):
    tot = arb(0)
    for m in range(MK+1):
        s = arb(2*m) + arb(1)/2
        tot += 2*((-s*a).exp() - (-s*b).exp())/s
    s0 = arb(2*(MK+1)) + arb(1)/2
    return tot + arb(0, 2*(-s0*a).exp()/(s0*(1 - (-2*a).exp())))
def H1s(a, b):
    tot = arb(0)
    for m in range(MK+1):
        s = arb(2*m) + arb(1)/2
        tot += 2*((-s*a).exp()*(s*a + 1) - (-s*b).exp()*(s*b + 1))/(s*s)
    s0 = arb(2*(MK+1)) + arb(1)/2
    return tot + arb(0, 2*(-s0*a).exp()*(a + 2/s0)/(s0*(1 - (-2*a).exp())))
EPS = arb('1e-6')
def _integ0(f, b, lead):
    if float(b.mid()) <= 1e-6: return lead(b)
    return lead(EPS) + acb.integral(f, acb(EPS), acb(b)).real
def R_int(b):
    def f(u, _): return u*(u/2).exp()/u.sinh() - 1 - u/2
    def lead(x): return -x**3/72 - x**4/64 + arb(0, 2*x**5)
    return _integ0(f, b, lead)
def R1_int(b):
    def f(u, _): return (u/2).exp()/u.sinh() - 1/u - arb(1)/2
    def lead(x): return -x**2/48 - x**3/48 + arb(0, 2*x**4)
    return _integ0(f, b, lead)
def H0f(a, b):
    if float(a.mid()) < 0.5: return (b/a).log() + (b - a)/2 + (R1_int(b) - R1_int(a))
    return H0s(a, b)
def H1f(a, b):
    if float(a.mid()) < 0.5: return (b - a) + (b*b - a*a)/4 + (R_int(b) - R_int(a))
    return H1s(a, b)
def H1_0(b): return b + b*b/4 + R_int(b)
_H03 = None; _W0 = None
def arch2(d):
    global _H03, _W0
    if _H03 is None:
        _H03 = sum((2*(-3*(arb(2*m)+arb(1)/2)).exp()/(arb(2*m)+arb(1)/2) for m in range(MK+1)), arb(0)) + arb(0, 1e-150)
        _W0 = (acb(arb(1)/4)).digamma().real - PI.log()
    b2d = beta2_val(arb(d))
    ad = abs(float(arb(d).mid()))
    pts = sorted(set([0.0, ad, 1-ad, 1+ad, 2-ad, 2+ad, 3.0]))
    pts = [p for p in pts if 0 <= p <= 3.0]
    tot = _W0*b2d
    for a, b in zip(pts[:-1], pts[1:]):
        if b - a < 1e-14: continue
        aa, bb = arb(str(a)), arb(str(b)); umid = (a + b)/2
        pp, qp = beta2_lin(arb(d), 1, umid)
        pm, qm = beta2_lin(arb(d), -1, umid)
        alpha = b2d - (pp + pm)/2; beta = -(qp + qm)/2
        if a == 0.0: tot += beta*H1_0(bb)
        else: tot += alpha*H0f(aa, bb) + beta*H1f(aa, bb)
    return tot + b2d*_H03
def prime2(d):
    dd = arb(d)
    return sum(((arb(p).log()/arb(n).sqrt())*(beta2_val(LOGN[n] + dd) + beta2_val(dd - LOGN[n])) for n, p in NS), arb(0))
def K2(d):
    dd = arb(d)
    return (BI2**2/2)*(dd/2).cosh() - prime2(dd) + arch2(dd)
# ---------------- shards -> A, H ----------------
def load_shard_dir():
    d = '/tmp/heart-src/experiments/round7_turbo'
    if glob.glob(d + '/shard_*.json'): return d
    out = '/tmp/flint15/shards'
    os.makedirs(out, exist_ok=True)
    data = subprocess.run(['git', '-C', os.path.expanduser('~/workspace/GPPDiscovery5'),
                           'archive', '1da974d', 'experiments/round7_turbo'],
                          capture_output=True, check=True).stdout
    tf = tarfile.open(fileobj=io.BytesIO(data))
    tf.extractall('/tmp/flint15/shardx')
    return '/tmp/flint15/shardx/experiments/round7_turbo'
def build_AH():
    sd = load_shard_dir()
    B = {}
    for fn in sorted(glob.glob(sd + '/shard_*.json')):
        for k, v in json.load(open(fn)).items():
            a, b = k.split(',')
            B[('d' if a == 'd' else int(a), 'd' if b == 'd' else int(b))] = arb(v)
    def Bph(j, k): return B[(j-1, k-1)] if j <= k else B[(k-1, j-1)]
    def Bpd(j): return B[(j-1, 'd')]
    g = B[('d', 'd')]
    basis = [('p', i) for i in range(1, 13)] + [('d', 0)] + [('s', k) for k in range(13, 73)]
    def Bval(x, y):
        def terms(z):
            if z[0] == 'p': return [(('p', z[1]), arb(1))]
            if z[0] == 'd': return [(('d', 0), arb(1))]
            return [(('p', z[1]), arb(1)), (('d', 0), -w_(z[1]))]
        tot = arb(0)
        for (zx, cx) in terms(x):
            for (zy, cy) in terms(y):
                if zx[0] == 'p' and zy[0] == 'p': v = Bph(zx[1], zy[1])
                elif zx[0] == 'p': v = Bpd(zx[1])
                elif zy[0] == 'p': v = Bpd(zy[1])
                else: v = g
                tot += cx*cy*v
        return tot
    A = [[Bval(basis[a], basis[b]) for b in range(73)] for a in range(73)]
    H = [[(4*PI*j)**2*(Bval(basis[a], ('p', j)) - w_(j)*Bval(basis[a], ('d', 0))) for j in range(1, 73)] for a in range(73)]
    return A, H, g
# ---------------- tail sums ----------------
def tail_sums():
    # phi1(j) = 1/om_j^2 ; phi2(j) = 1/((j-72) om_j^2) ; |w_j| = SQ8/om_j
    # sums over j>J0 by direct sum to J1 + integral bound for decreasing summands
    J1 = 200000
    S = {'w1': arb(0), 'w2': arb(0), '11': arb(0), '12': arb(0), '22': arb(0)}
    for j in range(J0+1, J1+1):
        wj = om(j); iw2 = 1/(wj*wj); ph2 = iw2/(j-72); aw = SQ8/wj
        S['w1'] += aw*iw2; S['w2'] += aw*ph2
        S['11'] += iw2*iw2; S['12'] += iw2*ph2; S['22'] += ph2*ph2
    # integral remainders: for decreasing f, sum_{j>J1} f(j) <= int_{J1}^inf f ; bound each by its J1-rate form:
    # w1 ~ c/j^3, w2 ~ c/j^4 (with (j-72) ~ j), 11 ~ c/j^4, 12 ~ c/j^5, 22 ~ c/j^6 ; use exact majorants:
    # |w|phi1 = SQ8/(16 pi^2 j^2 * 4 pi j) = SQ8/(64 pi^3 j^3): int = SQ8/(128 pi^3 J1^2)
    S['w1'] += SQ8/(128*PI**3*J1*J1) * arb('1.001')
    S['w2'] += SQ8/(64*PI**3) * (arb(1)/(3*(J1-72)**3)) * arb('1.34')   # (j/(j-72))^3 <= (200000/199928)^3 ~ 1.0011; 1.34 generous for lower terms
    S['11'] += arb(1)/(256*PI**4) * (arb(1)/(3*J1**3)) * arb('1.001')
    S['12'] += arb(1)/(256*PI**4) * (arb(1)/(4*(J1-72)**4)) * arb('1.34')
    S['22'] += arb(1)/(256*PI**4) * (arb(1)/(5*(J1-72)**5)) * arb('1.34')
    return S
def main():
    t00 = time.time()
    A, H, g = build_AH()
    print('A,H built', time.time()-t00, flush=True)
    q_delta = -(K2(arb(0)) - K2(arb(1)/2))/2
    print('q_delta =', q_delta, flush=True)
    # q vector on X' rows: a=0..11 phi_1..12 ; a=12 delta ; a=13..72 psi_k (k=a)
    q = [None]*73
    for i in range(1, 73):
        qphi = H[12][i-1]
        if i <= 12: q[i-1] = qphi
        if i >= 13: q[i] = qphi - w_(i)*q_delta
    q[12] = q_delta
    qabs = [abs(x.mid()) + x.rad() for x in q]
    # Q values for j=1..J0
    Qs = {}
    t0 = time.time()
    for j in range(1, J0+1):
        Qs[j] = Q_of(om(j))
    print(f'Q computed 1..{J0} in {time.time()-t0:.0f}s', flush=True)
    # C~ head blocks -> S_head
    S_head = arb_mat(73, 73, [arb(0)]*73*73)
    blk = 256
    for jstart in range(73, J0+1, blk):
        jend = min(jstart + blk, J0+1)
        cols = []
        for j in range(jstart, jend):
            wj = w_(j); Qj = Qs[j]
            rd = w4rho_delta(j, Qj)/om(j)**2
            crow = [None]*73
            cphi = {}
            for i in range(1, 73):
                cphi[i] = wj*H[12][i-1] + w4rho_phi(i, j, Qs[i], Qj)/om(j)**2
            for a in range(12): crow[a] = cphi[a+1]
            crow[12] = wj*q_delta + rd
            for k in range(13, 73): crow[k] = cphi[k] - w_(k)*(wj*q_delta + rd)
            cols.append(crow)
        M = arb_mat(73, jend-jstart, [cols[jj][a] for a in range(73) for jj in range(jend-jstart)])
        S_head = S_head + M*M.transpose()
        print(f'  block {jstart}..{jend-1} done {time.time()-t00:.0f}s', flush=True)
    # W_tail and rank-one tail
    Hsq = sum((arb(1)/(j*j) for j in range(1, J0+1)), arb(0))
    W_tail = (PI*PI/6 - Hsq)/(2*PI*PI)
    qqT = arb_mat(73, 73, [q[a]*q[b] for a in range(73) for b in range(73)])
    CCt = S_head + W_tail*qqT
    # tail entrywise radii
    S = tail_sums()
    Avec = [arb(0)]*73; Bvec = [arb(0)]*73
    for i in range(1, 73):
        if i <= 12: Bvec[i-1] = C_PHI
        else:
            Avec[i] = abs(w_(i))*C_DELTA; Bvec[i] = C_PHI
    Avec[12] = C_DELTA
    R_tail = [[None]*73 for _ in range(73)]
    maxrt = arb(0)
    for a in range(73):
        Ua = Avec[a]*S['w1'] + Bvec[a]*S['w2']
        for b in range(73):
            Ub = Avec[b]*S['w1'] + Bvec[b]*S['w2']
            Vab = Avec[a]*Avec[b]*S['11'] + (Avec[a]*Bvec[b] + Avec[b]*Bvec[a])*S['12'] + Bvec[a]*Bvec[b]*S['22']
            R_tail[a][b] = qabs[a]*Ub + qabs[b]*Ua + Vab
            if R_tail[a][b] > maxrt: maxrt = R_tail[a][b]
    print('max R_tail entry:', maxrt, flush=True)
    CCt = [[CCt[a, b] + arb(0, R_tail[a][b]) for b in range(73)] for a in range(73)]
    Hm = arb_mat(73, 72, [H[a][j] for a in range(73) for j in range(72)])
    HHt = Hm*Hm.transpose()
    Gram = [[HHt[a, b] + CCt[a][b] for b in range(73)] for a in range(73)]
    # ---- validations (float) ----
    import numpy as np
    Gm = np.array([[float(Gram[a][b].mid()) for b in range(73)] for a in range(73)])
    Grad = np.array([[float(Gram[a][b].rad()) for b in range(73)] for a in range(73)])
    Cm = np.array([[float(CCt[a][b].mid()) for b in range(73)] for a in range(73)])
    Crad = np.array([[float(CCt[a][b].rad()) for b in range(73)] for a in range(73)])
    Am = np.array([[float(A[a][b].mid()) for b in range(73)] for a in range(73)])
    res = {
      'q_delta': str(q_delta),
      'Gram_delta_delta': str(Gram[12][12]),
      'Gram_delta_delta_floatref_kgram': 0.0003352306686632827,
      'CCt_diag_rows_0_6_11_12_13_40_72': [float(Cm[i, i]) for i in [0, 6, 11, 12, 13, 40, 72]],
      'CCt_diag_ref_identity': [6.768918695580019e-08, 1.1176349821411328e-08, 3.438802161802912e-09, 9.805361317525695e-06, 1.0160422945570912e-10, 1.3511241315179982e-12, 2.606476406461366e-13],
      'max_abs_CCt': float(np.abs(Cm).max()),
      'Gram_max_radius': float(Grad.max()), 'CCt_max_radius': float(Crad.max()),
      'W_tail': str(W_tail),
    }
    Mmid = Am - Cm/5.24628
    Y = Mmid[12:, 12:]; E12 = Mmid[:12, :12]; C12 = Mmid[:12, 12:]
    Schur = E12 - C12 @ np.linalg.solve(Y, C12.T)
    res['float_schur_min_eig_lam5.24628'] = float(np.linalg.eigvalsh(Schur)[0])
    res['float_schur_ref_R10'] = 1.4116859583389672e-08
    print(json.dumps(res, indent=1), flush=True)
    # ---- dump balls (45-digit strings) ----
    def dump(M, fn):
        json.dump([[x.str(45) for x in row] for row in M], open(fn, 'w'))
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    dump(Gram, here + '/flint_gram_k_2026_11_15_balls.json')
    dump(CCt, here + '/flint_gram_k_2026_11_15_cct_balls.json')
    dump(A, here + '/flint_gram_k_2026_11_15_a_balls.json')
    json.dump(res, open(here + '/flint_gram_k_2026_11_15_results.json', 'w'), indent=1)
    print('TOTAL', time.time()-t00, flush=True)
if __name__ == '__main__':
    main()
