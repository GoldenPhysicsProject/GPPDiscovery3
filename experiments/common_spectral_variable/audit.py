"""Binary64 diagnostics, not interval certificates or a proof of RH."""
import cmath
import json
import math
from pathlib import Path


def psi(z):
    z = complex(z)
    shift = 0j
    while abs(z) < 32 or z.real < 16:
        shift -= 1 / z
        z += 1
    # Bernoulli asymptotic after recurrence into the right half-plane.
    coeff = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730]
    return shift + cmath.log(z) - 1/(2*z) - sum(
        b / (2*j*z**(2*j)) for j, b in enumerate(coeff, 1))


def gamma_m(k):
    return 1j*(psi(.25-.5j*k)-psi(.25))


def prime_m(p, k):
    q = cmath.exp(-(.5-1j*k)*math.log(p))
    return 1j*math.log(p)*(1+q)/(1-q)


def primes(n):
    a = bytearray(b'\1')*(n+1)
    a[:2] = b'\0\0'
    for p in range(2, math.isqrt(n)+1):
        if a[p]:
            a[p*p::p] = b'\0'*len(a[p*p::p])
    return [p for p in range(2, n+1) if a[p]]


def run():
    pts = [complex(x,y) for x in [-20,-3,0,3,20] for y in [.01,.25,1,4]]
    gm = min(gamma_m(k).imag for k in pts)
    pm = min(prime_m(p,k).imag for p in [2,3,5,101] for k in pts)
    assert gm > 0 and pm > 0
    residual = max(abs((prime_m(p,k)-1j*math.log(p))-
                       2j*math.log(p)/(cmath.exp((.5-1j*k)*math.log(p))-1))
                   for p in [2,3,5,101] for k in pts)
    # A local centered kernel already has a negative boundary value.
    centered_negative = (prime_m(2,math.pi/math.log(2))-1j*math.log(2)).imag
    assert centered_negative < 0
    orientation_residual = max(abs(prime_m(p,-k.conjugate())+
                                   prime_m(p,k).conjugate())
                               for p in [2,3,5,101] for k in pts)
    gamma_orientation_residual = max(abs(gamma_m(-k.conjugate())+
                                         gamma_m(k).conjugate()) for k in pts)
    ps = primes(100000)
    rows = []
    s = .75
    arch = (1/s+1/(s-1)-.5*math.log(math.pi)+.5*psi(s/2)).real
    for cutoff in [100,1000,10000,100000]:
        j = math.fsum(math.log(p)/(p**s-1) for p in ps if p <= cutoff)
        rows.append(dict(cutoff=cutoff, naive_completed_imag=arch-j,
                         first_order_scale=cutoff**(1-s)/(1-s)))
    result = dict(status='diagnostics only; exact proofs in private bridge note',
                  gamma_min_imag=gm, prime_min_imag=pm,
                  centering_identity_max_error=residual,
                  centered_p2_negative_boundary=centered_negative,
                  prime_orientation_error=orientation_residual,
                  gamma_orientation_error=gamma_orientation_residual,
                  cutoff_rows=rows)
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    run()
