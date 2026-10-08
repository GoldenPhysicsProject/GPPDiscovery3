"""Odd step-function Galerkin test of the completed Weil rank-one term.

No zeta zeros used. Binary64 quadrature, not a positivity certificate.
"""
import json
import math


def integral(f, a, b, tol=1e-11):
    def rec(a, b, fa, fm, fb, old, tol, depth):
        mid = (a+b)/2
        fl, fr = f((a+mid)/2), f((mid+b)/2)
        left = (mid-a)*(fa+4*fl+fm)/6
        right = (b-mid)*(fm+4*fr+fb)/6
        err = left+right-old
        if depth == 0 or abs(err) <= 15*tol:
            return left+right+err/15
        return rec(a, mid, fa, fl, fm, left, tol/2, depth-1) + rec(mid, b, fm, fr, fb, right, tol/2, depth-1)
    fa, fm, fb = f(a), f((a+b)/2), f(b)
    return rec(a, b, fa, fm, fb, (b-a)*(fa+4*fm+fb)/6, tol, 20)


def prime_powers(limit):
    sieve = bytearray(b'\x01')*(limit+1)
    sieve[:2] = b'\x00\x00'
    for p in range(2, math.isqrt(limit)+1):
        if sieve[p]:
            sieve[p*p:limit+1:p] = b'\x00'*((limit-p*p)//p+1)
    terms = []
    for p in range(2, limit+1):
        if sieve[p]:
            n = p
            while n <= limit:
                terms.append((math.log(n), math.log(p)/math.sqrt(n)))
                n *= p
    return terms


def cholesky(A):
    n = len(A)
    R = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            x = A[i][j]-math.fsum(R[i][k]*R[j][k] for k in range(j))
            if i == j:
                if x <= 0:
                    raise ValueError(('nonpositive pivot', i, x))
                R[i][j] = math.sqrt(x)
            else:
                R[i][j] = x/R[j][j]
    return R


def run(L, m):
    d = L/(2*m)
    rho = lambda y: math.exp(-y/2)/(-math.expm1(-2*y))
    # E(y) is linear between bin-shift knots, so integrate scalar weights.
    weights = [0.0]*(2*m+1)
    weights[1] = integral(lambda y: 1/(2*d) if y == 0 else y*rho(y)/d, 0, d)
    for k in range(1, 2*m):
        a, b = k*d, (k+1)*d
        weights[k] += integral(lambda y: rho(y)*(b-y)/d, a, b)
        weights[k+1] += integral(lambda y: rho(y)*(y-a)/d, a, b)
    primes = prime_powers(math.floor(math.exp(L)))
    for y, w in primes:
        pos = y/d
        k = min(int(pos), 2*m-1)
        u = pos-k
        weights[k] += w*(1-u)
        weights[k+1] += w*u
    # Unit-norm odd indicator pairs. Products of amplitudes are exactly +/-1/2.
    cells = [[(m+i, 1), (m-1-i, -1)] for i in range(m)]
    B = [[0.0]*m for _ in range(m)]
    for shift, weight in enumerate(weights):
        for i in range(m):
            for j in range(m):
                hp = math.fsum(si*sj/2 for ki, si in cells[i] for kj, sj in cells[j] if ki-kj == shift)
                hm = math.fsum(si*sj/2 for ki, si in cells[i] for kj, sj in cells[j] if kj-ki == shift)
                E = (2 if i == j else 0)-hp-hm
                B[i][j] += weight*E
    r = math.exp(-L/2)
    M_infty = math.log(8*math.pi)+0.5772156649015329+math.pi/2
    M_infty += math.log((1-r)/(1+r))-2*math.atan(r)
    M = M_infty+2*math.fsum(w for _, w in primes)
    for i in range(m):
        B[i][i] -= M
    v = [4/math.sqrt(2*d)*(math.cosh((i+1)*d/2)-math.cosh(i*d/2)) for i in range(m)]
    R = cholesky(B)
    z = []
    for i in range(m):
        z.append((v[i]-math.fsum(R[i][j]*z[j] for j in range(i)))/R[i][i])
    susceptibility = 2*math.fsum(x*x for x in z)
    Q = [[B[i][j]-2*v[i]*v[j] for j in range(m)] for i in range(m)]
    RQ = cholesky(Q)
    # Independently evaluate the original explicit formula on one odd vector.
    f = [(-1)**i/(i+1) for i in range(m)]
    N = math.fsum(x*x for x in f)
    def corr(y):
        pos = y/d
        k = min(int(pos), 2*m-1)
        u = pos-k
        def at(shift):
            return math.fsum(f[i]*f[j]*si*sj/2
                            for i in range(m) for j in range(m)
                            for ki, si in cells[i] for kj, sj in cells[j]
                            if ki-kj == shift)
        return (1-u)*at(k)+u*at(k+1)
    h1 = corr(d)
    def arch(y):
        if y == 0:
            return N/2+(h1-N)/d
        return (math.exp(y/2)*corr(y)-N)/math.sinh(y)
    arch_value = math.fsum(integral(arch, k*d, (k+1)*d, 1e-9) for k in range(2*m))
    arch_value += N*math.log(math.tanh(L/2))
    pole = -2*math.fsum(a*b for a, b in zip(v, f))**2
    direct = pole-(math.log(4*math.pi)+0.5772156649015329)*N-arch_value
    direct -= 2*math.fsum(w*corr(y) for y, w in primes)
    matrix = math.fsum(f[i]*Q[i][j]*f[j] for i in range(m) for j in range(m))
    assert abs(direct-matrix) < 2e-7
    return {'L': L, 'odd_bins': m, 'susceptibility': susceptibility,
            'margin_to_one': 1-susceptibility,
            'minimum_Q_cholesky_pivot': min(RQ[i][i]**2 for i in range(m)),
            'direct_vs_collision_error': abs(direct-matrix)}


if __name__ == '__main__':
    rows = [run(L, m) for L in (2, 4, 6, 8) for m in (8, 16)]
    print(json.dumps({'status': 'finite Galerkin binary64 checks, not an RH proof', 'rows': rows}, indent=2))
