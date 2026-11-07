#!/usr/bin/env python3
"""Extract exact generic nonzero-mu Ds=4 Yang-Mills sewing baselines.

Uses the convention-locked generic massive-vector tree engine.  The reconstruction
identity is C^(4)=C^(V_m)-C^(S), with vector and scalar sewings evaluated at the
same generic rational kinematics.  This script prints exact factorizations and
checks the threshold r=1 regression without making any box-only claim.

In physical variables beta=|p|/E, rho=mu/E and c=cos(theta), the same-helicity
channel collapses to

    C^(S)_same = rho^4/(1-beta*c)^2,
    Cv_same-3 Cs_same = 16 beta^2/(1-beta*c)^2,
    C^(4)_same = (2 rho^4 + 16 beta^2)/(1-beta*c)^2.

For the mixed-helicity channel set u=beta^2 sin(theta)^2.  Then

    C^(S)_mixed = u^2/(1-beta*c)^2,
    C^(4)_mixed = 2 (u^2-8u+8)/(1-beta*c)^2,
    C^(V_m)_mixed = (3u^2-16u+16)/(1-beta*c)^2.

These are exact symbolic identities in the rational chart, not fitted forms.
"""
from __future__ import annotations
import sympy as sp
import massive_vector_generic_state_sum_symbolic as g


def main() -> None:
    assert sp.simplify(g.beta**2 + g.rho**2 - 1) == 0

    Mmm = g.vector_tree_matrix(-1, -1)
    Mpp = g.vector_tree_matrix(+1, +1)
    Smm = g.scalar_tree(-1, -1)
    Spp = g.scalar_tree(+1, +1)
    Cv_same = g.sew(Mmm, Mpp)
    Cs_same = sp.factor(Smm * Spp)
    C4_same = sp.factor(sp.cancel(Cv_same - Cs_same))

    Mpm = g.vector_tree_matrix(+1, -1)
    Mmp = g.vector_tree_matrix(-1, +1)
    Spm = g.scalar_tree(+1, -1)
    Smp = g.scalar_tree(-1, +1)
    Cv_mixed = g.sew(Mpm, Mmp)
    Cs_mixed = sp.factor(Spm * Smp)
    C4_mixed = sp.factor(sp.cancel(Cv_mixed - Cs_mixed))

    same_defect = sp.factor(
        4 * (g.r**2 - 1)**2 * (g.t**2 + 1)**2 / (g.r**2 + g.t**2)**2
    )
    ctheta = (1 - g.t**2) / (1 + g.t**2)
    denom = (1 - g.beta * ctheta)**2

    physical_scalar_same = sp.factor(g.rho**4 / denom)
    physical_defect_same = sp.factor(16 * g.beta**2 / denom)
    physical_c4_same = sp.factor((2 * g.rho**4 + 16 * g.beta**2) / denom)
    beta_only_c4_same = sp.factor(
        2 * (g.beta**4 + 6 * g.beta**2 + 1) / denom
    )

    sin2 = sp.factor(1 - ctheta**2)
    u = sp.factor(g.beta**2 * sin2)
    physical_scalar_mixed = sp.factor(u**2 / denom)
    physical_c4_mixed = sp.factor(2 * (u**2 - 8*u + 8) / denom)
    physical_vector_mixed = sp.factor((3*u**2 - 16*u + 16) / denom)

    assert sp.simplify(C4_same - (Cv_same - Cs_same)) == 0
    assert sp.simplify(C4_mixed - (Cv_mixed - Cs_mixed)) == 0
    assert sp.simplify((Cv_same - 3 * Cs_same) - same_defect) == 0
    assert sp.simplify(C4_same - (2 * Cs_same + same_defect)) == 0

    assert sp.simplify(Cs_same - physical_scalar_same) == 0
    assert sp.simplify(same_defect - physical_defect_same) == 0
    assert sp.simplify(C4_same - physical_c4_same) == 0
    assert sp.simplify(C4_same - beta_only_c4_same) == 0

    assert sp.simplify(Cs_mixed - physical_scalar_mixed) == 0
    assert sp.simplify(C4_mixed - physical_c4_mixed) == 0
    assert sp.simplify(Cv_mixed - physical_vector_mixed) == 0

    assert sp.simplify((Cv_same - 3 * Cs_same).subs(g.r, 1)) == 0
    assert sp.simplify(same_defect.subs(g.r, 1)) == 0
    assert sp.simplify(Cv_mixed.subs(g.r, 1) - 16) == 0
    assert sp.simplify(Cs_mixed.subs(g.r, 1)) == 0
    assert sp.simplify(C4_mixed.subs(g.r, 1) - 16) == 0

    print("beta^2+rho^2 =", sp.simplify(g.beta**2 + g.rho**2))
    print("C4_same physical =", physical_c4_same)
    print("Cs_same physical =", physical_scalar_same)
    print("same defect physical =", physical_defect_same)
    print("u=beta^2 sin^2(theta) =", u)
    print("Cs_mixed physical =", physical_scalar_mixed)
    print("C4_mixed physical =", physical_c4_mixed)
    print("Cv_mixed physical =", physical_vector_mixed)
    print("threshold C4_same =", sp.factor(C4_same.subs(g.r, 1)))
    print("threshold C4_mixed =", sp.factor(C4_mixed.subs(g.r, 1)))
    print("PASS: exact generic Ds=4 physical-coordinate extraction")


if __name__ == "__main__":
    main()
