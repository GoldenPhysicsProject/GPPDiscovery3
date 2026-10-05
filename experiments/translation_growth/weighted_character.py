"""Numerical control for weighted continuity versus temperedness.

H_alpha = L2(R, exp(2 alpha |u|) du). The functional integrating exp(z u) f(u)
has exact dual norm sqrt(alpha/(alpha**2-delta**2)), delta=Re z, if |delta|<alpha.
This example contains no zeta zeros and establishes no arithmetic bound.
Run with Python 3; no external dependencies. Output is reproducible JSON.
"""
import cmath
import json
import math


def simpson(f, n):
    step = 2.0 / n
    return step / 3 * (f(-1.0) + f(1.0) + sum(
        (4 if k % 2 else 2) * f(-1.0 + k * step)
        for k in range(1, n)))


def bump(v):
    return math.exp(-1 / (1 - v*v)) if abs(v) < 1 else 0.0


def scan(n):
    alpha, delta, gamma = 0.5, 0.1, 3.0
    z = complex(delta, gamma)
    # b(v)=bump(v)*exp(-i gamma v) makes the base functional nonzero.
    base = simpson(lambda v: math.exp(delta*v)*bump(v), n)
    dual_norm = math.sqrt(alpha / (alpha*alpha - delta*delta))
    rows = []
    for t in [-80, -40, -10, 0, 10, 40, 80]:
        value = simpson(lambda v: cmath.exp(z*(t+v))*bump(v)*
                        cmath.exp(-1j*gamma*v), n)
        predicted = cmath.exp(z*t)*base
        norm = math.sqrt(simpson(lambda v: bump(v)**2 *
                                math.exp(2*alpha*abs(t+v)), n))
        error = abs(value-predicted)/abs(predicted)
        ratio = abs(value)/(dual_norm*norm)
        assert error < 1e-10, (t, error)
        assert ratio <= 1 + 1e-10, (t, ratio)
        rows.append(dict(t=t, amplitude=abs(value), weighted_norm=norm,
                         covariance_relative_error=error, cauchy_schwarz_ratio=ratio))
    return dict(alpha=alpha, delta=delta, gamma=gamma, panels=n,
                exact_dual_norm=dual_norm, base_value=base, rows=rows)


if __name__ == '__main__':
    coarse, fine = scan(1024), scan(2048)
    differences = [abs(a['amplitude']/b['amplitude']-1)
                   for a, b in zip(coarse['rows'], fine['rows'])]
    norm_differences = [abs(a['weighted_norm']/b['weighted_norm']-1)
                        for a, b in zip(coarse['rows'], fine['rows'])]
    assert max(differences+norm_differences) < 1e-8
    fine['resolution_relative_difference'] = max(differences+norm_differences)
    fine['status'] = 'floating-point control, not interval certification or an RH proof'
    print(json.dumps(fine, indent=2, sort_keys=True))
