#!/usr/bin/env python3
"""Exploratory full Wave234 S1 parity probe on exact-structure Arch annihilators.

The bump-train coefficients come from P_K(z)=prod_{m<K-1}(z-r_m),
r_m=exp(-0.7*(2*m+0.5)). All prime and pole formulas are the *source's*
Wave234 S1 normalization; the diagonal Arch integral is evaluated with
scipy only when installed. This is a floating-point diagnostic, NOT a
certified full-form sign. The Arch cross term uses its Laplace Gram sum.

Useful identities: the diagonal prime form uses autocorrelation
coefficients of P_K, while the sewing prime form uses coefficients of P_K^2.
Thus each is a short linear combination of compact-support scalar sums
T(h)=sum_{prime powers n} Lambda(n)n^-1/2 w_theta(n) G0(log n-h).
"""

import argparse
import bisect
import json
import math
from functools import lru_cache

A0 = 7 / 20
DELTA = 7 / 10
THETA = 1 / 4


def sine2_autocorrelation(x: float) -> float:
    """G0(x) = integral of b0(y)b0(y+x)dy, b0=sin^2(pi*y) on [0,1]."""
    t = abs(x)
    if t >= 1:
        return 0.0
    return ((1-t)*(2 + math.cos(2*math.pi*t)) / 8
            + 3*math.sin(2*math.pi*t)/(16*math.pi))


def sine2_laplace(z: float) -> float:
    """Integral_0^1 sin^2(pi*y)exp(z*y)dy (z real)."""
    if abs(z) < 1e-12:
        return 0.5
    pi2 = math.pi**2
    return 2*pi2*math.expm1(z)/(z*(z*z+4*pi2))


def annihilator(k: int) -> tuple[list[float], list[float]]:
    if k < 2:
        raise ValueError("K must be >=2")
    roots = [math.exp(-DELTA*(2*m+0.5)) for m in range(k-1)]
    coeff = [1.0]
    for root in roots:
        new = [0.0]*(len(coeff)+1)
        for j, val in enumerate(coeff):
            new[j] -= root*val
            new[j+1] += val
        coeff = new
    return roots, coeff


def polynomial_at_roots_product(x: float, roots: list[float]) -> float:
    return math.prod(x-r for r in roots)


@lru_cache(maxsize=4)
def prime_power_atoms(limit: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Prime powers <=limit as sorted (log(n), Lambda(n)/sqrt(n)) atoms."""
    sieve = bytearray(b"\x01")*(limit+1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, math.isqrt(limit)+1):
        if sieve[p]:
            count = (limit-p*p)//p+1
            sieve[p*p:limit+1:p] = b"\x00"*count
    atoms = []
    for p in range(2, limit+1):
        if not sieve[p]:
            continue
        lp = math.log(p)
        power = p
        while power <= limit:
            atoms.append((math.log(power), lp/math.sqrt(power)))
            if power > limit//p:
                break
            power *= p
    atoms.sort()
    return tuple(t[0] for t in atoms), tuple(t[1] for t in atoms)


def prime_kernel(logs: tuple[float,...], weights: tuple[float,...],
                 shift: float, theta: float) -> float:
    start = bisect.bisect_left(logs, shift-1)
    stop = bisect.bisect_right(logs, shift+1)
    if theta == 0:
        return math.fsum(weights[i]*sine2_autocorrelation(logs[i]-shift)
                         for i in range(start,stop))
    return math.fsum(weights[i]*(2*math.cosh(theta*logs[i]))
                     *sine2_autocorrelation(logs[i]-shift)
                     for i in range(start,stop))


def prime_parts(coeff: list[float], theta: float) -> dict:
    k = len(coeff)
    # Largest sewing support: log(n) <=2(a0+delta*(K-1))+2.
    maxlog = 2*(A0 + DELTA*(k-1)) + 2
    cutoff = math.ceil(math.exp(maxlog))
    logs, weights = prime_power_atoms(cutoff)
    # D = sum_d corr_d*T(-delta*d), C=sum_l conv_l*T(2a0+delta*l+1).
    corr = [math.fsum(coeff[j]*coeff[j+d] for j in range(k-d))
            for d in range(k)]
    conv = [math.fsum(coeff[j]*coeff[l-j]
             for j in range(max(0,l-k+1), min(k-1,l)+1))
            for l in range(2*k-1)]
    kernel = lambda h: prime_kernel(logs, weights, h, theta)
    diag = corr[0]*kernel(0)
    # Positive d comes from v_j v_{j+d} twice (±d), but T(+d) and T(-d)
    # differ because support log(n)>=0.
    diag += math.fsum(corr[d]*(kernel(DELTA*d)+kernel(-DELTA*d))
                     for d in range(1,k))
    cross = math.fsum(conv[l]*kernel(2*A0+DELTA*l+1)
                      for l in range(2*k-1))
    return {"2D": 2*diag, "X": cross, "atom_count":len(logs),
            "prime_cutoff":cutoff}


def pole_parts(roots: list[float], theta: float) -> dict:
    poles = [0.0, 1.0] if theta == 0 else [1+theta, 1-theta, theta, -theta]
    plus = []
    minus = []
    for p in poles:
        z = p-0.5
        aplus = math.exp(z*A0)*sine2_laplace(z)*polynomial_at_roots_product(
            math.exp(DELTA*z), roots)
        aminus = math.exp(-z*A0)*sine2_laplace(-z)*polynomial_at_roots_product(
            math.exp(-DELTA*z), roots)
        plus.append((aplus+aminus)**2)
        minus.append((aplus-aminus)**2)
    return {"plus":math.fsum(plus)/2, "minus":-math.fsum(minus)/2,
            "poles":poles}


def continuum_sewing_prime(roots: list[float], theta: float) -> float:
    """Exact continuum-main-term identity evaluated in floats.

    The atomic weight Lambda(n)/sqrt(n) has density exp(t/2) dt under the
    prime number theorem. The F_theta weight adds 2*cosh(theta*t).
    Since C_V(t)=(V*V)(t), the continuum substitution factorizes into
    squares of the corresponding one-sided Mellin/Laplace amplitudes.
    This is a convolution identity, NOT a bound for discrete prime errors.
    """
    shifts = [0.0] if theta == 0 else [-theta,theta]
    return math.fsum((math.exp((.5+s)*A0)*sine2_laplace(.5+s)
                      *polynomial_at_roots_product(math.exp(DELTA*(.5+s)),roots))**2
                     for s in shifts)


def gaussian_q_binomial_edge(k: int) -> list[float]:
    """Coefficients of the reversed P_K via the finite q-binomial theorem.

    Return [v_{K-1}, v_{K-2}, ..., v_0]. q=exp(-2*DELTA), r0=exp(-DELTA/2).
    The j-th coefficient is (-r0)^j q^{j(j-1)/2} [K-1 choose j]_q.
    """
    m = k-1
    q = math.exp(-2*DELTA)
    r0 = math.exp(-DELTA/2)
    result = []
    for j in range(k):
        qb = math.prod((1-q**(m-l))/(1-q**(l+1)) for l in range(j))
        result.append((-r0)**j*q**(j*(j-1)//2)*qb)
    return result


def folded_cross_arch(roots: list[float], theta: float,
                      max_modes: int = 80) -> float:
    """Return v*C_arch*v <=0. Positive in source's opposite-Gram sign."""
    shifts = [0.0] if theta == 0 else [-theta,theta]
    vals = []
    for n in range(max_modes):
        for shift in shifts:
            lam = 2*n+0.5+shift
            moment = (math.exp(-A0*lam)*sine2_laplace(-lam)
                      *polynomial_at_roots_product(math.exp(-DELTA*lam),roots))
            vals.append(moment*moment)
    return -math.fsum(vals)


def diagonal_arch(roots: list[float], theta: float, cutoff: float) -> tuple[float,float]:
    """Source gamma/digamma Fourier integral, numerical scipy quad only."""
    from scipy.integrate import quad
    from scipy.special import digamma

    def integrand(t):
        if abs(t) < 1e-10:
            sinc = 0.5
        elif abs(t-2*math.pi) < 1e-8:
            sinc = -0.25
        else:
            sinc = 4*math.pi**2*math.sin(t/2)/(t*(4*math.pi**2-t*t))
        prod = math.prod(1+r*r-2*r*math.cos(DELTA*t) for r in roots)
        if theta == 0:
            m = digamma(0.25+0.5j*t).real-math.log(math.pi)
        else:
            m = (digamma(0.25-theta/2+0.5j*t).real
                 +digamma(0.25+theta/2+0.5j*t).real-2*math.log(math.pi))
        return m*sinc*sinc*prod/math.pi

    # Separate removable t=2pi and sign-change region from oscillatory tail.
    segments = [0.0, 2*math.pi, 20.0, 50.0, 100.0, 200.0, cutoff]
    segments = sorted(set(t for t in segments if t<=cutoff))
    terms, errors = [], []
    for lo,hi in zip(segments,segments[1:]):
        val,err = quad(integrand,lo,hi,epsabs=8e-13,epsrel=2e-12,limit=600)
        terms.append(val)
        errors.append(err)
    return math.fsum(terms), math.fsum(errors)


def sector_matrix_check(k: int, theta: float, cutoff: float) -> dict:
    """Independent 2 parity matrices from Q Toeplitz and C Hankel formulas.

    This independently checks all source normalizations against Turbo's
    Wave234 eigenvalue table. Numpy/scipy are optional exploration tools.
    """
    import numpy as np
    from scipy.integrate import quad
    from scipy.special import digamma

    maxlog = 2*(A0+DELTA*(k-1))+2
    logs, weights = prime_power_atoms(math.ceil(math.exp(maxlog)))
    kernel = lambda h:prime_kernel(logs,weights,h,theta)
    qprime = [kernel(DELTA*d)+kernel(-DELTA*d) for d in range(k)]
    cprime = [kernel(2*A0+DELTA*l+1) for l in range(2*k-1)]

    def arch_measure(t):
        if abs(t) < 1e-10:
            envelope = .5
        elif abs(t-2*math.pi) < 1e-8:
            envelope = -.25
        else:
            envelope = 4*math.pi**2*math.sin(t/2)/(t*(4*math.pi**2-t*t))
        if theta == 0:
            m = digamma(.25+.5j*t).real-math.log(math.pi)
        else:
            m = (digamma(.25-theta/2+.5j*t).real
                 +digamma(.25+theta/2+.5j*t).real-2*math.log(math.pi))
        return m*envelope**2/math.pi

    arch_toeplitz = []
    cuts = sorted(set(t for t in
                      (0,2*math.pi,20,50,100,200,cutoff) if t<=cutoff))
    for d in range(k):
        arch_toeplitz.append(math.fsum(
            quad(lambda t:arch_measure(t)*math.cos(DELTA*d*t),lo,hi,
                 epsabs=8e-13,epsrel=2e-12,limit=600)[0]
            for lo,hi in zip(cuts,cuts[1:])))

    shifted = [0.0] if theta == 0 else [-theta,theta]
    arch_hankel = [0.0]*(2*k-1)
    for l in range(2*k-1):
        s = 2*A0+DELTA*l
        arch_hankel[l] = -math.fsum(
            math.exp(-(2*n+.5+shift)*s)*
            sine2_laplace(-(2*n+.5+shift))**2
            for n in range(80) for shift in shifted)

    poles = [0.,1.] if theta==0 else [1+theta,1-theta,theta,-theta]
    plus, minus = np.zeros((k,k)),np.zeros((k,k))
    for p in poles:
        z=p-.5
        xp=np.array([math.exp(z*(A0+DELTA*j))*sine2_laplace(z)
                     for j in range(k)])
        xm=np.array([math.exp(-z*(A0+DELTA*j))*sine2_laplace(-z)
                     for j in range(k)])
        plus += .5*np.outer(xp+xm,xp+xm)
        minus -= .5*np.outer(xp-xm,xp-xm)

    for j in range(k):
        for h in range(k):
            diagonal_arch_prime=arch_toeplitz[abs(j-h)]-qprime[abs(j-h)]
            cross_arch_prime=arch_hankel[j+h]-cprime[j+h]
            plus[j,h]+=diagonal_arch_prime+cross_arch_prime
            minus[j,h]+=diagonal_arch_prime-cross_arch_prime
    _, coeff = annihilator(k)
    vector=np.array(coeff)
    denom=float(vector@vector)
    return {"lambda_min_plus":float(np.linalg.eigvalsh(plus)[0]),
            "lambda_min_minus":float(np.linalg.eigvalsh(minus)[0]),
            "annihilator_rayleigh_plus":float((vector@plus@vector)/denom),
            "annihilator_rayleigh_minus":float((vector@minus@vector)/denom),
            "max_symmetry_residual":float(max(np.max(abs(plus-plus.T)),
                                              np.max(abs(minus-minus.T))))}


def evaluate(k:int, theta:float, arch_cutoff:float|None) -> dict:
    roots, coeff = annihilator(k)
    norm = math.fsum(c*c for c in coeff)
    pole = pole_parts(roots,theta)
    prime = prime_parts(coeff,theta)
    cross = folded_cross_arch(roots,theta)
    continuum=continuum_sewing_prime(roots,theta)
    prime_residual=prime["X"]-continuum
    report = {"K":k,"theta":theta,"norm2":norm,"prime":prime,
              "pole":pole,"cross_arch":cross,"vector":coeff,
              "cross_arch_over_norm":cross/norm,
              "continuum_prime_sewing":continuum,
              "arithmetic_prime_residual":prime_residual,
              "gaussian_q_binomial_error":max(abs(v-w) for v,w in zip(
                  reversed(coeff),gaussian_q_binomial_edge(k)))}
    if theta==0:
        report["exact_pole_continuum_identity_error"]=(pole["plus"]-continuum)
    if arch_cutoff is not None:
        diag, quadrature_err = diagonal_arch(roots,theta,arch_cutoff)
        report["diagonal_arch"] = diag
        report["quadrature_error_reported"] = quadrature_err
        report["arch_t_cutoff"] = arch_cutoff
        report["parity"] = {}
        for sign,label in [(+1,"plus"),(-1,"minus")]:
            p = pole[label]
            a = diag+sign*cross
            n = prime["2D"]+sign*prime["X"]
            report["parity"][label]={
                "pole_over_norm":p/norm,
                "arch_over_norm":a/norm,
                "prime_over_norm":n/norm,
                "full_over_norm":math.fsum([p,a,-n])/norm}
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k",type=int,default=6)
    ap.add_argument("--theta",type=float,default=0,
                    choices=(0,0.25))
    ap.add_argument("--arch-cutoff",type=float,default=None,
                    help="include scipy diagonal-Arch integral through t cutoff")
    ap.add_argument("--matrix-check",action="store_true",
                    help="also assemble all parity matrices for Turbo source audit")
    args = ap.parse_args()
    if args.k>10:
        ap.error("K<=10 protects single-workstation prime sieve memory")
    output=evaluate(args.k,args.theta,args.arch_cutoff)
    if args.matrix_check:
        if args.arch_cutoff is None:
            ap.error("--matrix-check requires --arch-cutoff")
        output["matrix_check"]=sector_matrix_check(args.k,args.theta,args.arch_cutoff)
        for label in ("plus","minus"):
            p=output["parity"][label]["full_over_norm"]
            q=output["matrix_check"][f"annihilator_rayleigh_{label}"]
            output["matrix_check"][f"projection_residual_{label}"]=p-q
    print(json.dumps({
        "status":"NUMERICAL diagnostics only; gamma integral truncated; no RH proof",
        "source":"GPP-bridge/MuseResearch/turbo-source-wave234-s1-schur-definitions.md",
        "output":output},indent=2))


if __name__ == "__main__":
    main()
