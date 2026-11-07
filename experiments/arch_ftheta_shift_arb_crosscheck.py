#!/usr/bin/env python3
"""Independent python-flint Arb review of the Ftheta shifted Arch witness.

This evaluates the actual one shifted Gram square at 384-bit precision,
independently of the integer fixed-point interval implementation.
The mathematical proof still requires its exact Laplace/Gram identities;
Arb is an additional independent arithmetic implementation check.
"""

from flint import arb, ctx, fmpq

from arch_cross_superexp_null import exact_upper
from arch_ftheta_shift_gap import theta_mode_rayleigh_lower


def as_arb(q):
    return arb(fmpq(q.numerator, q.denominator))


def result(k: int) -> None:
    a0 = arb(fmpq(7, 20))
    delta = arb(fmpq(7, 10))
    lam_star = arb(fmpq(1, 4))
    pi2 = arb.pi() ** 2

    def s(lam):
        return 2*pi2*(1-(-lam).exp()) / (lam*(lam*lam+4*pi2))

    roots = [(-delta*arb(fmpq(4*n+1, 2))).exp() for n in range(k-1)]
    coeff = [arb(1)]
    for r in roots:
        p = [arb(0)] * (len(coeff)+1)
        for j, c in enumerate(coeff):
            p[j] -= c*r
            p[j+1] += c
        coeff = p
    norm = sum((c*c for c in coeff), arb(0))
    assert norm > 0
    shifted_root = (-delta*lam_star).exp()
    value = arb(1)
    for r in roots:
        value *= shifted_root-r
    shifted_one_mode = ((-a0*lam_star).exp()*s(lam_star)*value)**2/norm
    lower = theta_mode_rayleigh_lower(k)
    assert shifted_one_mode > as_arb(lower), (k, str(shifted_one_mode), str(lower))

    # An independent finite evaluation checks that the undisturbed head
    # is consistent with, but does NOT itself prove, the infinite tail bound.
    zeta_head = arb(0)
    for n in range(k-1, k+20):
        lam = arb(fmpq(4*n+1, 2))
        z = (-delta*lam).exp()
        pval = arb(1)
        for r in roots:
            pval *= z-r
        term = (-a0*lam).exp()*s(lam)*pval
        zeta_head += term*term
    zeta_head /= norm
    assert zeta_head < as_arb(exact_upper(k))
    print(f"K={k}: 384-bit Arb one-mode theta lower and zeta head agree "
          "with independent rational certificates; PASS")


def main() -> None:
    ctx.prec = 384
    for k in (6, 10):
        result(k)


if __name__ == "__main__":
    main()
