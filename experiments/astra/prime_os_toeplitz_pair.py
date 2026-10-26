"""Prime-field OS/KMS Fredholm phase test, including reciprocal orientations.

This identifies a real operator-theoretic transition at Re(s)=1/2; it
does NOT assert that zeta zeros coincide with Toeplitz eigenmodes.

The half-density local character is z_p(s)=p**(1/2-s), and the canonical
unilateral shift S on l2(N0) gives D_p(s)=S-z_p(s)I. Its adjoint has the
twisted TFD/KMS occupation vector in its kernel iff Re(s)>1/2.

This is a standard Toeplitz/Fredholm calculation applied to arithmetic
half-density characters and the opposite-orientation reciprocal pair.
"""

from __future__ import annotations

import cmath
import math

import numpy as np


def one_orientation(p: int, sigma: float, t: float, n: int = 80) -> dict:
    s=complex(sigma,t)
    z=cmath.exp((0.5-s)*math.log(p))
    modulus=abs(z)
    # Rectangular finite-section embedding (n+1 x n): unlike square
    # compressions it retains the last shifted coordinate, so it preserves
    # the true infinite-space lower bound ||(S-z)x||>=(1-|z|)||x||.
    # Rectangular matrices DO NOT reveal the infinite Fredholm index;
    # index below follows from the analytic shift theorem.
    d=np.zeros((n+1,n),dtype=np.complex128)
    for j in range(n):
        d[j,j]=-z
        d[j+1,j]=1
    singular=np.linalg.svd(d,compute_uv=False)
    # Infinite-space statements, from ||Sx||=||x|| and Toeplitz index.
    index=-1 if modulus<1-1e-12 else (0 if modulus>1+1e-12 else None)
    bound=abs(1-modulus)
    # For |z|<1, the exact cokernel vector of D is proportional to
    # (1,conj(z),conj(z)^2,...) and is normalized by sqrt(1-|z|**2).
    return {"p":p,"sigma":sigma,"t":t,"abs_z":modulus,
            "fredholm_index":index,"infinite_lower_gap_bound":bound,
            "finite_min_singular":float(min(singular)),
            "local_os_vacuum_normalizable":modulus<1-1e-12}


def paired_orientation(p: int,sigma: float,t: float,n: int=80):
    # s and 1-s, not complex conjugation. The physical on-axis pairing
    # satisfies z_p(1-s)=1/z_p(s), so at sigma=1/2 this is conj(z_p(s)).
    plus=one_orientation(p,sigma,t,n)
    minus=one_orientation(p,1-sigma,-t,n)
    assert math.isclose(plus["abs_z"]*minus["abs_z"],1,rel_tol=1e-12)
    return plus,minus


if __name__=="__main__":
    first_zeta_height=14.134725141734693
    for sigma in (0.3,0.45,0.5,0.55,0.7):
        x,y=paired_orientation(2,sigma,first_zeta_height)
        print(f"sigma={sigma:.2f} plus |z|={x['abs_z']:.8f}, "
              f"index={x['fredholm_index']}, min_sv80={x['finite_min_singular']:.7g}; "
              f"shadow |z|={y['abs_z']:.8f}, index={y['fredholm_index']}, "
              f"min_sv80={y['finite_min_singular']:.7g}")
        assert sigma==0.5 or x["fredholm_index"]!=y["fredholm_index"]
        if sigma==0.5:
            assert x["fredholm_index"] is None and y["fredholm_index"] is None
