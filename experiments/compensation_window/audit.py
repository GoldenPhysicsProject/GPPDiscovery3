"""Exploratory checks of fixed-window radial compensation; not an RH proof.

Standard-library binary64 arithmetic, with an analytic series-tail bound.
No interval certificate or Lean verification is claimed.
"""
import bisect
import json
import math


def main():
    ell = math.log(2)
    xmax = 100000
    limit = 2 * xmax
    sieve = bytearray(b'\x01') * (limit + 1)
    sieve[0:2] = b'\x00\x00'
    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p:limit+1:p] = b'\x00' * ((limit-p*p)//p+1)
    powers = {}
    for p in range(2, limit+1):
        if sieve[p]:
            n = p
            while n <= limit:
                powers[n] = math.log(p)
                n *= p
    ns = sorted(powers)
    logs = [math.log(n) for n in ns]
    weights = [powers[n]/math.sqrt(n) for n in ns]

    def prefixes(ws):
        w, v = [0.0], [0.0]
        for z, a in zip(logs, ws):
            w.append(w[-1]+a)
            v.append(v[-1]+a*z)
        return w, v

    ordinary = prefixes(weights)
    correction = prefixes([w/n for n, w in zip(ns, weights)])

    def triangle(x, pref):
        t = math.log(x)
        lo, mid, hi = (bisect.bisect_right(ns, y) for y in (x/2, x, 2*x))
        w, v = pref
        return ((ell-t)*(w[mid]-w[lo]) + v[mid]-v[lo]
                +(ell+t)*(w[hi]-w[mid])-v[hi]+v[mid])/ell

    A = (6*math.sqrt(2)-8)/ell
    knots = {2.0, float(xmax)}
    for n in ns:
        knots.update(x for x in (n/2, float(n), 2.0*n) if 2 <= x <= xmax)
    largest = (-math.inf, None)
    for x in sorted(knots):
        deficit = A*math.sqrt(x)-triangle(x, ordinary)
        if deficit > largest[0]:
            largest = (deficit, x)
    samples = []
    for x in (2, 10, 100, 1000, 10000, 100000):
        r = triangle(x, correction)
        bound = math.exp(3*ell/2)*(math.exp(ell)-math.exp(-ell)+1/x)
        bound *= (math.log(x)+ell)/math.sqrt(x)
        direct = math.fsum(w/n*max(0, 1-abs(math.log(n/x))/ell)
                           for n, w in zip(ns, weights) if x/2 <= n <= 2*x)
        assert abs(r-direct) < 1e-9
        assert -1e-9 <= r <= bound
        samples.append({'x': x, 'R': r, 'elementary_upper_bound': bound})
    N = 100000
    terms = []
    for k in range(N):
        a, b = 2*k+.5, 2*k+1
        terms.append(.5/(a*b)+math.expm1(-a*ell)/(ell*a*a))
    c0_partial = 2*A-math.log(4*math.pi)-0.5772156649015329-2*math.fsum(terms)
    tail = (.5+1/ell)/(2*(N-1))
    print(json.dumps({
        'arithmetic': 'binary64; no rigorous rounding certificate',
        'xmax': xmax, 'knots_tested': len(knots),
        'largest_sampled_deficit': {'value': largest[0], 'x': largest[1]},
        'correction_samples': samples,
        'C0_partial': c0_partial, 'analytic_tail_bound': tail,
        'C0_plus_Dbound_estimated_upper': c0_partial+tail+math.pi**2/(24*ell)
    }, indent=2))


if __name__ == '__main__':
    main()
