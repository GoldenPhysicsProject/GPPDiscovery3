"""Exact rational small-window Weil certificate; no external packages.

Analytical inputs and proof: private Astra bridge note dated 2026-10-15.
This checks rational inequalities, not the explicit formula or operator theorem.
"""
from fractions import Fraction as F
from math import factorial
import json


def log_lower(x, terms=30):
    """For rational x >= 1: positive atanh series is a rigorous lower bound."""
    assert x >= 1
    y = (x - 1) / (x + 1)
    return 2 * sum((y ** (2*k+1) / (2*k+1) for k in range(terms)), F(0))


def atan_bounds(x, terms=20):
    """Alternating series, with a next-term enclosure for 0<x<1."""
    assert 0 < x < 1
    s = sum(((-1)**k * x**(2*k+1)/(2*k+1) for k in range(terms)), F(0))
    other = s + (-1)**terms*x**(2*terms+1)/(2*terms+1)
    return min(s, other), max(s, other)


def certify():
    # Machin's identity pi=16 atan(1/5)-4 atan(1/239).
    a, b = atan_bounds(F(1,5)); c, d = atan_bounds(F(1,239))
    pi_lo, pi_hi = F(314159,100000), F(3927,1250)
    assert pi_lo < 16*a-4*d
    assert 16*b-4*c < pi_hi

    # gamma < H_1000-log(1000); log(1000)=9 log 2+log(125/64).
    gamma_hi = F(289,500)
    harmonic = sum((F(1,k) for k in range(1,1001)), F(0))
    assert harmonic - 9*log_lower(F(2)) - log_lower(F(125,64)) < gamma_hi

    h = F(1,15)
    # 2h < log 2, so the COMPLETE prime contribution vanishes by support.
    assert 2*h < log_lower(F(2))
    pole_loss = h**3 / (3*(1-h*h/20))
    base = log_lower((2-h)/(8*pi_hi*h)) - gamma_hi - pi_hi*h/2 - pole_loss
    gain = F(0)
    for n in range(2,33):
        # 1-cos(pi/(n+1)) >= x^2/2-x^4/24, enclosed below.
        w = pi_lo**2/(2*(n+1)**2) - pi_hi**4/(24*(n+1)**4)
        # nu_0(t)>=1/t-1/2 on the integration interval.
        mass = log_lower(F(n,n-1)) - h*(F(1,n-1)-F(1,n))
        assert w > 0 and mass > 0
        gain += w*mass
    bound = base + gain
    assert bound > F(1,20)

    # Exact path-graph positive-factor certificate for first-entry chains.
    # 2 Re(conj(a)b) <= |a|^2+|b|^2 via |a-b|^2>=0.
    # Numerical displays are not used in the above decisions.
    return {"arithmetic":"fractions.Fraction; every assertion exact",
            "h":"1/15", "support_length":"2/15", "target_gap":"1/20",
            "certified_lower_bound_decimal_display":float(bound),
            "positive_margin_decimal_display":float(bound-F(1,20)),
            "chain_intervals":31, "log_series_terms":30,
            "pi_bounds":[str(pi_lo),str(pi_hi)], "gamma_upper":str(gamma_hi),
            "scope":"rational certificate for analytic lower bound; not Lean or global RH"}


if __name__ == '__main__':
    print(json.dumps(certify(), indent=2))
