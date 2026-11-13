#!/usr/bin/env python3
"""FS3 corollary probe (oracle advice 2026-11-13, section C; ORACLE SIGHTING,
unaudited — this probe is the audit of its checkable corollary).

FS3 claims: D(c) = 0 iff M_c == 0, where M_c(s) = sum_i c_i (e^{(s-1/2)x_i} +
e^{-(s-1/2)x_i}); on the line M_c(1/2 + it) = 2 sum_i c_i cos(t x_i).
Checkable corollary: the joint-minimisation near-null configurations
(smallest banked: lam_min = +6.06e-9, extension run) should have M_c ~ 0
measurably.

Probe: at the extension's Nelder-Mead minimiser, take the extremal c of the
defect form (proj_eigmin), and measure sup_t and RMS_t of |M_c(1/2+it)|,
normalised by 2||c||_2, over t in [0,500] and [0,2000]. Compare with the
same statistics at a generic configuration (equispaced L=2, n=4). For
distinct frequencies the RMS ratio tends to 1/sqrt(2) = 0.7071 by
almost-periodicity, whatever c is; a near-annihilator configuration must
beat that decisively. Also report the coefficient pattern (near-duplicate
point groups) so the mechanism of any smallness is visible.

NUMERICAL. Writes flint_fs3_probe_2026_11_13_results.json.
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import research.flint_archdom_scout_2026_11_12 as scout  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def mc_stats(x, c, tmax_list):
    x = np.asarray(x, float)
    c = np.asarray(c, float)
    out = {}
    for tmax in tmax_list:
        t = np.linspace(0.0, tmax, int(tmax / 0.001) + 1)
        Mc = 2.0 * (np.cos(np.outer(t, x)) @ c)
        nrm = 2.0 * float(np.linalg.norm(c))
        out[str(tmax)] = {
            "sup_ratio": float(np.max(np.abs(Mc)) / nrm),
            "rms_ratio": float(np.sqrt(np.mean(Mc ** 2)) / nrm),
        }
    return out


def main():
    ext_res = json.load(open(os.path.join(
        HERE, "flint_archdom_margin_extension_2026_11_12_results.json")))
    rows = ext_res["minimisation_n4"]
    best = min(rows, key=lambda r: r["lam_min"])
    x_min = best["x"]
    L = best["L"]
    A, PmS, Pm, D, K, md = scout.build_mats(np.array(x_min), L)
    lam, lam_max, c = scout.proj_eigmin(D, np.array(x_min))
    print(f"NM minimiser: L={L} lam_min(banked)={best['lam_min']:.6e} "
          f"recomputed={lam:.6e} spread={best['spread']:.4f}")
    print(f"  x = {x_min}")
    print(f"  extremal c = {[round(float(v), 6) for v in c]}  (||c||=1)")
    print(f"  ell_c(c) = {float(np.cosh(np.array(x_min)/2) @ c):.3e} (ker check)")

    xg = scout.configs(2.0, 4)["equispaced"]
    Ag, PmSg, Pmg, Dg, Kg, mdg = scout.build_mats(xg, 2.0)
    lamg, _, cg = scout.proj_eigmin(Dg, xg)
    print(f"\ngeneric (equispaced L=2, n=4): lam_min={lamg:.6e}")
    print(f"  extremal c = {[round(float(v), 6) for v in cg]}")

    res = {
        "minimiser": {"x": x_min, "L": L, "lam_min_banked": best["lam_min"],
                      "lam_min_recomputed": float(lam),
                      "c_extremal": [float(v) for v in c],
                      "mc_stats": mc_stats(x_min, c, [500.0, 2000.0])},
        "generic": {"x": [float(v) for v in xg], "lam_min": float(lamg),
                    "c_extremal": [float(v) for v in cg],
                    "mc_stats": mc_stats(xg, cg, [500.0, 2000.0])},
        "reference": "RMS ratio for distinct frequencies -> 1/sqrt(2) = 0.7071",
    }
    print("\nM_c statistics (ratios to 2||c||_2):")
    for name in ["minimiser", "generic"]:
        for tmax, st in res[name]["mc_stats"].items():
            print(f"  {name:>9} t<={tmax}: sup={st['sup_ratio']:.4f} rms={st['rms_ratio']:.4f}")
    ms = res["minimiser"]["mc_stats"]["500.0"]
    gs = res["generic"]["mc_stats"]["500.0"]
    res["verdict"] = (
        "FS3 corollary SUPPORTED in uniform form on [0,500] (NUMERICAL): the "
        "near-null minimiser's M_c is measurably small where the generic "
        "configuration's is O(1)"
        if ms["sup_ratio"] < 0.35 * gs["sup_ratio"] else
        "FS3 corollary FAILS in uniform form (NUMERICAL): the near-null "
        "minimiser's M_c is NOT uniformly small; its smallness lives only in "
        "the pairing (D itself), which restates the phenomenon instead of "
        "explaining it")
    print(f"\n{res['verdict']}")
    with open(os.path.join(HERE, "flint_fs3_probe_2026_11_13_results.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("wrote flint_fs3_probe_2026_11_13_results.json")


if __name__ == "__main__":
    main()
