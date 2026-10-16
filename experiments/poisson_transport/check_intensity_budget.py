"""Floating-point diagnostic for Astra's finite Poisson transport.
Not an interval certificate and not evidence of global Weil positivity/RH.
The exact derivation is recorded in Astra's private research notes.
"""
import json
import math

def prime_powers(limit):
    prime = [True] * (limit + 1)
    prime[:2] = [False, False]
    for p in range(2, limit + 1):
        if not prime[p]:
            continue
        for multiple in range(p*p, limit + 1, p):
            prime[multiple] = False
        n = p
        while n <= limit:
            yield n, math.log(p)
            n *= p

def diagnostic(limit=100):
    intensity = 0.0
    rows = []
    for n, lam in sorted(prime_powers(limit)):
        u = math.log(n)
        w = lam / math.sqrt(n)
        alpha = 2 * lam / (n - 1)
        intensity += alpha
        mean = alpha * math.sinh(u / 2)
        entropy = u * w - 2 * alpha * (math.cosh(u / 2) - 1)
        rows.append(dict(n=n, reference_intensity=intensity,
                         mean_residual=mean-w,
                         entropy=entropy, entropy_budget=u*w,
                         fits_reference=intensity <= 4))
    return rows

if __name__ == "__main__":
    print(json.dumps(diagnostic(), indent=2))
