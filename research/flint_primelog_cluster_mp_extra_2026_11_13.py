#!/usr/bin/env python3
"""
Supplementary mpmath dps=40 confirmation for flint_primelog_cluster
(Flint / Discovery 3, 2026-11-13). NUMERICAL.

The main run confirmed at dps40 its five smallest float cases (all Family B,
all positive). This supplement confirms the remaining small-margin cases the
main run did not reach:
  * Family A, eps=0.05, n=6 at the bracketing centres P=101 and P=1000003
    (float lam_min +4.138e-6, the family's smallest);
  * Family B at P=1009, n=4/6/8 (float lam_min +4e-9 / ~0 / ~0);
  * Family C at P=1009 (float lam_min +4.1403e-5).
Same projected-spectrum construction as the main run's mp section.
"""
import json
import math
import sys
import time

import numpy as np

sys.path.insert(0, ".")
import research.flint_archdom_scout_2026_11_12 as scout  # noqa: E402
import research.flint_archdom_margin_extension_2026_11_12 as ext  # noqa: E402
from research.flint_primelog_cluster_2026_11_13 import first_primes_from  # noqa: E402

import mpmath as mp
mp.mp.dps = 40


def triG_mp(y):
    ay = abs(y)
    if ay <= 1:
        return mp.mpf(2) / 3 - y * y + ay ** 3 / 2
    if ay <= 2:
        return (2 - ay) ** 3 / 6
    return mp.mpf(0)


def prime_mp(d):
    d = mp.mpf(d)
    lim = mp.e ** (2 + abs(d))
    s = mp.mpf(0)
    for n, lp in ext.PP:
        if n > lim:
            break
        ln = mp.log(n)
        s += mp.mpf(lp) / mp.sqrt(n) * (triG_mp(ln + d) + triG_mp(-ln + d))
    return s


def arch_mp(d):
    d = mp.mpf(d)
    f = lambda t: (mp.mpf(1) if t == 0 else (mp.sin(t / 2) / (t / 2)) ** 2) ** 2 \
        * mp.cos(t * d) * (mp.re(mp.digamma(mp.mpc(mp.mpf("0.25"), t / 2))) - mp.log(mp.pi))
    if d >= 2:
        nper = int(100 * d / mp.pi) + 1
        pts = [mp.mpf(0)] + [mp.mpf(k) * mp.pi / d for k in range(1, nper + 1)] + [mp.inf]
        return (1 / mp.pi) * mp.quad(f, pts)
    return (1 / mp.pi) * mp.quad(f, [0, mp.inf])


_amp, _pmp = {}, {}
def Amp(d):
    k = mp.nstr(mp.mpf(d), 15)
    if k not in _amp:
        _amp[k] = arch_mp(mp.mpf(d))
    return _amp[k]
def Pmp(d):
    k = mp.nstr(mp.mpf(d), 15)
    if k not in _pmp:
        _pmp[k] = prime_mp(mp.mpf(d))
    return _pmp[k]


bI2 = 16 * mp.sinh(mp.mpf("0.25")) ** 2


def proj_lam_min_mp(xlist):
    xs = [mp.mpf(v) for v in xlist]
    n = len(xs)
    ell = [mp.cosh(v / 2) for v in xs]
    sv = [mp.sinh(v / 2) for v in xs]
    Dm = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            d = abs(xs[i] - xs[j])
            Dm[i, j] = Amp(d) - Pmp(d) - 2 * bI2 ** 2 * sv[i] * sv[j]
    cols = []
    for k in range(1, n):
        v = mp.matrix(n, 1)
        v[k] = 1
        v[0] = -ell[k] / ell[0]
        cols.append(v)
    Z = []
    for v in cols:
        w = mp.matrix(v)
        for q in Z:
            w -= q * (q.T * w)[0]
        Z.append(w / mp.sqrt((w.T * w)[0]))
    m = n - 1
    Mp_ = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            Mp_[i, j] = (Z[i].T * Dm * Z[j])[0]
    ev = mp.eigsy(Mp_, eigvals_only=True)
    return min(ev)


if __name__ == "__main__":
    t0 = time.time()
    cases = []
    for P in [101, 1000003]:
        u0 = math.log(P)
        x = [u0 + v for v in np.linspace(-0.05, 0.05, 6)]
        cases.append((f"A P={P} eps=0.05 n=6", x))
    for n in [4, 6, 8]:
        ps = first_primes_from(1009, n)
        cases.append((f"B P=1009 n={n} primes={ps[0]}..{ps[-1]}",
                      [math.log(p) for p in ps]))
    offs = [-math.log(3), -math.log(2), 0.0, math.log(2), math.log(3)]
    cases.append(("C P=1009 offsets=log(1/3..3)", [math.log(1009) + v for v in offs]))
    out = {"label": "NUMERICAL", "experiment": "flint_primelog_cluster_mp_extra", "rows": []}
    for name, x in cases:
        mn = proj_lam_min_mp(x)
        print(f"{name}: mp lam_min = {mp.nstr(mn, 15)}")
        out["rows"].append({"case": name, "lam_min_mp": mp.nstr(mn, 20),
                            "verdict": "SURVIVES (NUMERICAL)" if mn >= 0 else "NEGATIVE at dps40 (NUMERICAL)"})
    out["runtime_s"] = time.time() - t0
    with open("research/flint_primelog_cluster_mp_extra_2026_11_13_results.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"done ({out['runtime_s']:.1f}s)")
