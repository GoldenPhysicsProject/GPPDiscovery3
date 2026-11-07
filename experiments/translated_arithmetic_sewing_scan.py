#!/usr/bin/env python3
"""Complete source-side translated prime--Arch interference for one K window.

This differs from the previously sampled *zero-side* translation control:
every point is independently computed from the finite, exactly supported
prime-power sum, the PNT main term, the Archimedean Laplace modes, and the
fixed diagonal gamma integral.  The exact identities for V_a(u)=V(u-a) are

    E(a) = Arch_diag(V) - 2 D_prime(V) = E(0),
    X(a) = sum_n Lambda(n)/sqrt(n) (V*V)(log n - 2a),
    Pole_plus(a) = exp(a) A^2, Pole_minus(a) = -exp(a) A^2,
    C_arch(a) = -sum_m exp(-2 lambda_m a) M_m(V)^2,
    P_plus(a) = E + C_arch(a) - [X(a)-exp(a) A^2],
    P_minus(a) = E - C_arch(a) + [X(a)-exp(a) A^2].

The gamma integral and floating cancellation are NUMERICAL, not rigorous.
Finite first-N known critical-line zeros can optionally cross-check C(a);
neither the scan nor those known zeros establish RH or an all-a bound.
"""

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arch_annihilator_full_s1_probe import (
    A0, DELTA, annihilator, continuum_sewing_prime, diagonal_arch,
    folded_cross_arch, prime_kernel, prime_parts, prime_power_atoms,
    sine2_laplace,
)


def convolution_coefficients(coeff: list[float]) -> list[float]:
    k = len(coeff)
    return [math.fsum(coeff[j]*coeff[l-j]
                      for j in range(max(0, l-k+1), min(k-1, l)+1))
            for l in range(2*k-1)]


def translated_prime_cross(logs, weights, conv, a: float) -> float:
    """Finite exact-support prime-power sum for sewn X(V_a), in floats."""
    return math.fsum(c*prime_kernel(logs, weights,
                                    2*(A0+a)+DELTA*l+1, 0.)
                     for l,c in enumerate(conv))


def translated_cross_arch(roots: list[float], a: float,
                          max_modes: int = 100) -> float:
    """Full folded cross Arch as an exponentially converging mode sum."""
    return -math.fsum(
        (math.exp(-(A0+a)*lam)*sine2_laplace(-lam)
         *math.prod(math.exp(-DELTA*lam)-r for r in roots))**2
        for lam in (2*m+.5 for m in range(max_modes)))


def known_zero_interference(k: int, zeros: int, shifts: list[float]) -> tuple:
    """Only first `zeros` known line pairs; not an all-zero identity check."""
    import mpmath as mp
    mp.mp.dps = 45
    a0,delta=mp.mpf(7)/20,mp.mpf(7)/10
    roots=[mp.exp(-delta*(2*m+mp.mpf('.5'))) for m in range(k-1)]
    base=[]
    for n in range(1,zeros+1):
        z=1j*mp.im(mp.zetazero(n))
        s=2*mp.pi**2*mp.expm1(z)/(z*(z*z+4*mp.pi**2))
        f=mp.exp(a0*z)*s*mp.fprod(mp.exp(delta*z)-r for r in roots)
        base.append((z, f*f))
    interference=[float(2*mp.re(mp.fsum(coef*mp.exp(2*a*z)
                                      for z,coef in base))) for a in shifts]
    # On-line known pairs have diagonal weights 2|F(i gamma)|^2.
    # This table is NOT an assertion that every zeta zero is on the line.
    modes=[{"n":j+1,"gamma":float(mp.im(z)),
            "diagonal_weight":float(2*abs(coef))}
           for j,(z,coef) in enumerate(base)]
    return interference,modes


def run(k: int, shifts: list[float], arch_cutoff: float,
        max_prime_cutoff: int, zeros: int = 0) -> dict:
    if k < 2 or not shifts or min(shifts) < 0:
        raise ValueError("K>=2 and a nonempty nonnegative translation grid required")
    roots,coeff=annihilator(k)
    # Support of (V_a*V_a): [2(A0+a), 2(A0+a+delta*(k-1)+1)].
    prime_cutoff=math.ceil(math.exp(2*(A0+max(shifts)+DELTA*(k-1)+1)))
    if prime_cutoff > max_prime_cutoff:
        raise ValueError(f"Prime sieve needs {prime_cutoff:,}; cap={max_prime_cutoff:,}")
    logs,weights=prime_power_atoms(prime_cutoff)
    conv=convolution_coefficients(coeff)
    continuum=continuum_sewing_prime(roots, 0.)
    # E is translation-invariant; only this one fixed quadrature is needed.
    arch,err=diagonal_arch(roots,0.,arch_cutoff)
    qprime=prime_parts(coeff,0.)["2D"]
    even=arch-qprime
    spectral,modes=known_zero_interference(k,zeros,shifts) if zeros else (None,[])
    rows=[]
    for i,a in enumerate(shifts):
        x=translated_prime_cross(logs,weights,conv,a)
        pole=math.exp(a)*continuum
        residual=math.fsum([x,-pole])
        carch=translated_cross_arch(roots,a)
        sewing=math.fsum([carch,-residual])
        plus=math.fsum([even,sewing])
        minus=math.fsum([even,-sewing])
        row={"a":a,"prime_X":x,"pnt_pole":pole,
             "epsilon_X_minus_pole":residual,"cross_arch":carch,
             "even_E":even,"sewing_C":sewing,"plus":plus,"minus":minus,
             "relative_saturation":abs(sewing)/abs(even) if even else None}
        if spectral is not None:
            row["known_line_zeros_C"]=spectral[i]
            row["known_zeros_source_residual"]=sewing-spectral[i]
        rows.append(row)
    a0_cross=translated_cross_arch(roots,0.)
    cross_regression=abs(a0_cross-folded_cross_arch(roots,0.))
    result={"K":k,"translations":shifts,"atom_count":len(logs),
            "prime_cutoff":prime_cutoff,"even_E":even,
            "quadrature_reported_error":err,"arch_cutoff":arch_cutoff,
            "fixed_cross_arch_regression":cross_regression,
            "min_parity_gap":min(min(row["plus"],row["minus"]) for row in rows),
            "max_relative_saturation":max(row["relative_saturation"] for row in rows),
            "zeros_pairs_optional":zeros,"rows":rows}
    if modes:
        result["known_zero_first_10_modes"]=modes[:10]
        result["known_zero_first_mode_fraction_of_even"]=modes[0]["diagonal_weight"]/even
        result["known_zero_first_two_fraction_of_even"]=sum(
            mode["diagonal_weight"] for mode in modes[:2])/even
        result["known_zeros_max_source_difference"]=max(
            abs(row["known_zeros_source_residual"]) for row in rows)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--k",type=int,default=6)
    p.add_argument("--a-max",type=float,default=2.)
    p.add_argument("--step",type=float,default=.2)
    p.add_argument("--arch-cutoff",type=float,default=500.)
    p.add_argument("--max-prime-cutoff",type=int,default=8_000_000)
    p.add_argument("--zeros",type=int,default=0,
                   help="optional first known-line zero pairs (mpmath)")
    args=p.parse_args()
    if args.step<=0 or args.a_max<0 or args.k>10:
        p.error("step>0, a-max>=0, 2<=K<=10 required")
    n=int(args.a_max/args.step+1e-12)
    shifts=[round(j*args.step,12) for j in range(n+1)]
    if args.zeros and not 10<=args.zeros<=400:
        p.error("optional known-line zero count must be 10..400")
    result=run(args.k,shifts,args.arch_cutoff,
               args.max_prime_cutoff,args.zeros)
    if args.zeros>=120:
        assert result["known_zeros_max_source_difference"] < 2e-8, (
            "Independent finite zero/source numerical check failed",
            result["known_zeros_max_source_difference"])
    print(json.dumps({"status":"NUMERICAL full source scan; not a uniform bound or RH proof",
                      "result":result},indent=2))


if __name__=="__main__":
    main()
