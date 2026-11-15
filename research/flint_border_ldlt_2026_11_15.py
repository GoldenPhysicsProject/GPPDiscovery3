# Flint run 15, item 4 (2026-11-15): 600-bit ball LDL^T of M = A_X' - CC^T/lambda_Z - beta*P
# in R7's scaled basis (E_12 by 1/sqrt(1e-8), delta by 1/sqrt(g), psi_k by D_k^{-1/2},
# D_k = 3 w(om_k)/om_k^4), elimination order [psi_13..72, delta, E_12].
# Certificate (R7's diag(r) argument): with M = M0 + E, |E_ab| <= R_ab entrywise and
# r_a = sum_b R_ab (scaled), LDL^T of M0 - diag(r) - beta*P~ completing with every pivot
# enclosed > 0 implies Schur_{E_12}(A_X' - CC^T/lambda_Z) >= beta for every matrix in the
# input balls. Inputs: flint_gram_k_2026_11_15_a_balls.json / _cct_balls.json (item 3').
# lambda_Z = 5.24628 (LOCKED, Heart (T1')); primary target beta2 = 1.3e-8 (pre-named);
# strengthened probe at 1.40e-8 also certified; 1.41e-8 fails (true float Schur 1.41171e-8).
from flint import arb, acb, ctx
import json, time, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
LAM = arb('5.24628')
def load(fn):
    return [[arb(s) for s in row] for row in json.load(open(fn))]
A = load(here + '/flint_gram_k_2026_11_15_a_balls.json')
C = load(here + '/flint_gram_k_2026_11_15_cct_balls.json')
Am = np.array([[float(A[a][b].mid()) for b in range(73)] for a in range(73)])
Cm = np.array([[float(C[a][b].mid()) for b in range(73)] for a in range(73)])
Mf = Am - Cm/float(LAM.mid())
Y = Mf[12:, 12:]; E12 = Mf[:12, :12]; C12 = Mf[:12, 12:]
Schur = E12 - C12 @ np.linalg.solve(Y, C12.T)
float_schur = float(np.linalg.eigvalsh(Schur)[0])
print('float Schur min eig:', float_schur, flush=True)
ctx.prec = 256
PI = arb.pi()
g_mid = A[12][12].mid()
s = [None]*73
for a in range(12): s[a] = arb(1)/10000
s[12] = g_mid.sqrt()
for k in range(13, 73):
    w = 4*PI*k
    what = (acb(arb(1)/4) + acb(0,1)*(w/2)).digamma().real - PI.log()
    Dk = 3*what.mid()/(w.mid()**4)
    s[k] = Dk.sqrt()
ctx.prec = 600
M = [[A[a][b] - C[a][b]/LAM for b in range(73)] for a in range(73)]
Mt = [[M[a][b]/(s[a]*s[b]) for b in range(73)] for a in range(73)]
r = [sum((Mt[a][b].rad() for b in range(73)), arb(0)) for a in range(73)]
perm = list(range(13, 73)) + [12] + list(range(12))
n = 73
def certify(beta):
    N = [[Mt[perm[i]][perm[j]].mid() for j in range(n)] for i in range(n)]
    for i in range(n):
        N[i][i] = N[i][i] - r[perm[i]]
        if perm[i] < 12:
            N[i][i] = N[i][i] - beta/(s[perm[i]]*s[perm[i]])
    L = [[arb(0)]*n for _ in range(n)]
    d = [arb(0)]*n
    for k in range(n):
        dk = N[k][k] - sum((L[k][i]*L[k][i]*d[i] for i in range(k)), arb(0))
        d[k] = dk
        if not (float(dk.lower()) > 0):
            return False, d, k
        for j in range(k+1, n):
            L[j][k] = (N[j][k] - sum((L[j][i]*L[k][i]*d[i] for i in range(k)), arb(0)))/dk
    return True, d, None
res = {'lambda_Z': '5.24628', 'float_schur_min_eig': float_schur,
       'max_scaled_radius_rowsum': max(float(x.mid()) for x in r), 'levels': {}}
for name, beta in [('beta2_prenamed', arb('1.3e-8')), ('strengthened', arb('1.40e-8')), ('probe_fail', arb('1.41e-8'))]:
    t0 = time.time()
    ok, d, badk = certify(beta)
    res['levels'][name] = {'beta': str(beta), 'certified': ok,
                           'min_pivot_lower': (min(float(x.lower()) for x in d) if ok else None),
                           'E_block_pivots_lower': ([float(d[k].lower()) for k in range(61, 73)] if ok else None),
                           'failed_pivot': badk, 'elapsed_s': time.time()-t0}
    if name == 'beta2_prenamed':
        res['pivots_beta2'] = [str(x) for x in d]
    print(name, beta, '->', 'CERTIFIED' if ok else f'FAILED at pivot {badk}', flush=True)
json.dump(res, open(here + '/flint_border_ldlt_2026_11_15_results.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != 'pivots_beta2'}, indent=1))
