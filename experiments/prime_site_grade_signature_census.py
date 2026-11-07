#!/usr/bin/env python3
"""Per-prime / formal-grade signature census in Wave234 S1 normalization.

For P^sigma = Pole^sigma + Toeplitz_arch - Toeplitz_prime
              + sigma*(Hankel_arch - Hankel_prime), isolate EXACTLY its
signed contribution from a SINGLE finite prime p. No pole or gamma
contribution is assigned to any individual finite prime.

The analytic expressions are the source's sine² compact autocorrelation
and positive prime-power weights; eigen-inertia decisions are FLOATING
and require an independent directed interval gate to become proofs.

A *single prime-power* changes by a positive scalar under Ftheta, so its
signature does not change. A full prime SITE combines many powers,
which theta weights unequally: its matrix signature can change.
"""
import argparse
import json
import math

import sys
from pathlib import Path

# Direct script invocation by rh_lab.runner puts experiments/ on sys.path;
# pytest/unittest module invocation starts at the repository root.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from arch_annihilator_full_s1_probe import (
    A0, DELTA, sine2_autocorrelation,
)


def site_atoms(p: int, k: int) -> list[tuple[int, float, float]]:
    """Return (exponent, log(n), Lambda(n)/sqrt(n)), n=p^exponent."""
    if p < 2 or k < 2:
        raise ValueError("p>=2 and K>=2")
    maxlog = 2*(A0 + DELTA*(k-1)) + 2
    logp = math.log(p)
    return [(a, a*logp, logp*math.exp(-a*logp/2))
            for a in range(1, math.ceil(maxlog/logp)+1)
            if a*logp < maxlog]


def site_matrix(p: int, k: int, sector: int, theta: float,
                *, include_conjugate: bool = False) -> list[list[float]]:
    """Finite-prime contribution to P^sector, no pole/arch pieces.

    theta=0 -> ordinary zeta weight, one grade.
    theta=+/-1/4 -> either shifted factor separately.
    include_conjugate -> *sum* of +/-theta shifted factors.
    """
    if sector not in (-1, 1):
        raise ValueError("sector must be +/-1")
    atoms = site_atoms(p, k)

    def kernel(h: float) -> float:
        return math.fsum(
            w * (2*math.cosh(theta*t) if include_conjugate
                 else math.exp(theta*t))
            * sine2_autocorrelation(t-h)
            for _,t,w in atoms
        )
    d = [kernel(DELTA*j)+kernel(-DELTA*j) for j in range(k)]
    c = [kernel(2*A0+DELTA*l+1) for l in range(2*k-1)]
    return [[-(d[abs(i-j)]+sector*c[i+j]) for j in range(k)]
            for i in range(k)]


def leading_minors_3(matrix: list[list[float]]) -> list[float]:
    if len(matrix) != 3 or any(len(row) != 3 for row in matrix):
        raise ValueError("currently supported exact inertia formula: K=3")
    a,b,c=matrix[0]
    d,e,f=matrix[1]
    g,h,i=matrix[2]
    return [a,a*e-b*d,a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)]


def inertia_from_ldl_3(matrix: list[list[float]]) -> dict:
    """A 3x3 symmetric matrix with nonzero leading minors is congruent
    to diag(Delta1, Delta2/Delta1, Delta3/Delta2). Floating diagnostic.
    """
    ds=leading_minors_3(matrix)
    assert max(abs(matrix[i][j]-matrix[j][i])
               for i in range(3) for j in range(3)) < 1e-12
    pivots=[ds[0],ds[1]/ds[0],ds[2]/ds[1]]
    if min(abs(x) for x in pivots) < 1e-9:
        return {"status":"too_close_to_zero", "minors":ds,"pivots":pivots}
    return {"status":"numerical", "positive":sum(x>0 for x in pivots),
            "negative":sum(x<0 for x in pivots),"minors":ds,"pivots":pivots}


def run(p: int=3, sector: int=-1) -> dict:
    k=3
    mats={
        "zeta_weight_0":site_matrix(p,k,sector,0),
        "shift_weight_plus_half":site_matrix(p,k,sector,.25),
        "shift_weight_minus_half":site_matrix(p,k,sector,-.25),
        "Ftheta_combined":site_matrix(
            p,k,sector,.25,include_conjugate=True),
    }
    return {
        "status":"FLOATING finite prime-site signature diagnostic; not a proof",
        "prime":p,"K":k,"parity":sector,
        "source":"Wave234 S1 signed prime blocks, pole and Arch omitted",
        "num_powers":len(site_atoms(p,k)),
        "formal_weights_are_not_geometric_grades":True,
        "matrices":{key:{"matrix":m,**inertia_from_ldl_3(m)}
                    for key,m in mats.items()},
        "warning":"Do not call any per-place block a complete centred form; "
                  "pole and Arch are global and the Ftheta grading is formal."
    }


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prime",type=int,default=3)
    ap.add_argument("--sector",type=int,choices=[-1,1],default=-1)
    args=ap.parse_args()
    print(json.dumps(run(args.prime,args.sector),indent=2))
