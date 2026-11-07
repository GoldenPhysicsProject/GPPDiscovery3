#!/usr/bin/env python3
"""Numerical audit of the 4D MHV Yang--Mills and gravity cuts.

All external momenta are outgoing. The cut routing is
  p1+p2+l1+l2=0,
  -l1-l2+p3+p4=0.

Checks, in the stripped conventions of YM_GRAVITY_MHV_TWO_PARTICLE_CUT.md,

  C_YM/A_tree = -i s t/(D1 D2),

and, using four-point KLT tree by tree,

  C_GR/M_tree = i s^3 t u/(D1 D2 D3 D4).

Several center-of-mass energies are included so the s^3 gravity scaling is
actually tested. This is an audit, not the analytic proof.
"""
import numpy as np


def pmat(p):
    E, px, py, pz = p
    return np.array([[E + pz, px - 1j * py],
                     [px + 1j * py, E - pz]], dtype=complex)


def spinors(p):
    M = pmat(p)
    i, j = np.unravel_index(np.argmax(np.abs(M)), M.shape)
    lam = M[:, j].copy()
    if abs(lam[i]) < 1e-13:
        raise RuntimeError("degenerate factorization pivot")
    tlam = np.zeros(2, dtype=complex)
    tlam[j] = 1.0
    tlam[1 - j] = M[i, 1 - j] / lam[i]
    if np.max(np.abs(np.outer(lam, tlam) - M)) > 1e-10:
        raise RuntimeError("spinor factorization failed")
    return lam, tlam


def angle(a, b):
    return a[0] * b[1] - a[1] * b[0]


def parke_taylor(lams):
    num = angle(lams[0], lams[1]) ** 4
    den = 1.0 + 0j
    for i in range(4):
        den *= angle(lams[i], lams[(i + 1) % 4])
    return 1j * num / den


def mdot(a, b):
    return a[0] * b[0] - np.dot(a[1:], b[1:])


def msq(a):
    return mdot(a, a)


def amp(order):
    return parke_taylor([spinors(p)[0] for p in order])


def gravity_tree(order):
    # M4(1,2,3,4) = -i s12 A(1,2,3,4) A(1,2,4,3)
    s12 = msq(order[0] + order[1])
    return -1j * s12 * amp(order) * amp([order[0], order[1], order[3], order[2]])


def one_case(E, theta, phi):
    p1 = np.array([-E, 0.0, 0.0, -E])
    p2 = np.array([-E, 0.0, 0.0,  E])
    p3 = np.array([ E, E*np.sin(theta), 0.0, E*np.cos(theta)])
    p4 = np.array([ E,-E*np.sin(theta), 0.0,-E*np.cos(theta)])
    l1 = np.array([ E, E*np.sin(phi), 0.0, E*np.cos(phi)])
    l2 = np.array([ E,-E*np.sin(phi), 0.0,-E*np.cos(phi)])

    assert np.max(np.abs(p1+p2+p3+p4)) < 1e-12
    assert np.max(np.abs(p1+p2+l1+l2)) < 1e-12
    assert np.max(np.abs(-l1-l2+p3+p4)) < 1e-12

    # Yang--Mills first color ordering.
    Aext = amp([p1, p2, p3, p4])
    AL = amp([p1, p2, l2, l1])
    AR = amp([-l1, -l2, p3, p4])
    ym_lhs = AL * AR / Aext

    s = msq(p1+p2)
    t = msq(p2+p3)
    u = msq(p1+p3)
    D1 = msq(l1+p1)
    D2 = msq(l1-p4)
    D3 = msq(l1+p2)
    D4 = msq(l1-p3)
    ym_rhs = -1j * s * t / (D1 * D2)

    # Gravity from KLT on each of the three four-point trees.
    Mext = gravity_tree([p1,p2,p3,p4])
    ML = gravity_tree([p1,p2,l2,l1])
    MR = gravity_tree([-l1,-l2,p3,p4])
    gr_lhs = ML * MR / Mext
    gr_rhs = 1j * s**3 * t * u / (D1 * D2 * D3 * D4)
    return ym_lhs, ym_rhs, gr_lhs, gr_rhs


def main():
    worst_ym = 0.0
    worst_gr = 0.0
    for E in (0.7, 1.0, 1.3):
        for theta in (0.5, 0.9, 1.4):
            for phi in (0.3, 0.7, 1.2, 2.0):
                if abs(theta-phi) < 0.05:
                    continue
                yl, yr, gl, gr = one_case(E, theta, phi)
                ey = abs(yl-yr) / max(1.0, abs(yr))
                eg = abs(gl-gr) / max(1.0, abs(gr))
                worst_ym = max(worst_ym, ey)
                worst_gr = max(worst_gr, eg)
                assert ey < 3e-12, (E, theta, phi, yl, yr, ey)
                assert eg < 5e-11, (E, theta, phi, gl, gr, eg)
    print("PASS YM: C_s/A_tree = -i s t/(D1 D2)")
    print("PASS GR: C_s/M_tree = i s^3 t u/(D1 D2 D3 D4)")
    print("worst YM relative error:", worst_ym)
    print("worst GR relative error:", worst_gr)


if __name__ == "__main__":
    main()
