"""
flint_inputb_generators_audit_2026_11_14.py -- second-arm audit of the input (b)
minorant generator chain (GPP-bridge work/claude-code/rung1: b_exact_mirror.py
+ cert_data.pieces_G + gen_lean_minorant.py) against the Lean pieceTable on
GPPVerify PR #597 head c930c5e (GppWindowMinorantL036).

Chain being audited, per piece i (centre c_i = (2i+1)/8, u = r - c_i, |u| <= 1/8):
  (1) Lean pieceTable[i]  ==  r_coeffs(pieces_G()[i])            [exact fmpq]
  (2) p_i(u) <= lb_i(u) - err_i  on [-1/8, 1/8]                  [exact fmpq,
      Bernstein-basis certificate: all Bernstein coefficients of
      (lb - err - p) on the interval are >= 0]
  (3) lb_i(u) - err_i <= G_80(c_i + u) at dense samples          [Arb balls,
      prec 200; the analytic Taylor-remainder step, checked numerically --
      its Lean proof is the remaining input-(b) content]
Steps (1)+(2) are exact; step (3) is NUMERICAL evidence about the mirror's
remainder inequalities (no mis-signing, no mis-scaling), not a proof.

Flint, 2026-11-14. Run: ~/workspace/venvs/arb-venv/bin/python (python-flint),
with the bridge rung1 directory on the path for the generator imports.
"""

import json
import re
import sys
from math import comb

sys.path.insert(0, "/home/hatch/workspace/gpp-bridge-flintwork/work/claude-code/rung1")

import flint
from flint import arb, ctx, fmpq

ctx.prec = 200

import b_exact_mirror as bm
import cert_data as cd

# ---------- (1) Lean table == regenerated G-based minorant ----------
src = open("/tmp/WindowMinorantL036.lean").read()
m = re.search(r"def pieceTable : List \(List ℚ\) := \[(.*?)\n\]", src, re.S)
fracs = re.findall(r"\((-?\d+) / (\d+) : ℚ\)", m.group(1))
lean_pieces = [[fmpq(int(n), int(d)) for (n, d) in fracs[13 * i:13 * (i + 1)]]
               for i in range(100)]

gen = cd.pieces_G()
mismatch = []
for (i, c, p_u, _eps) in gen:
    p_r = cd.r_coeffs(c, p_u, fmpq(0))
    if p_r != lean_pieces[i]:
        mismatch.append(i)
print(f"(1) Lean pieceTable == pieces_G regeneration: "
      f"{'EXACT MATCH on all 100 pieces' if not mismatch else 'MISMATCH on pieces ' + str(mismatch)}")


# ---------- helpers ----------
def to_u_coeffs_from_r(p_r, c):
    """ascending r-coefficients -> ascending (r - c) coefficients, exact."""
    out = [fmpq(0)] * len(p_r)
    for k, a in enumerate(p_r):
        for j in range(k + 1):
            out[j] += a * comb(k, j) * c ** (k - j)
    return out


def bernstein_min(coeffs_u, h):
    """Min Bernstein coefficient of poly(u) (ascending u-coeffs) on [-h, h].
    Substitute u = -h + 2h t, t in [0,1]; if all Bernstein coeffs >= 0 the
    polynomial is >= 0 on the interval."""
    N = len(coeffs_u) - 1
    a = [fmpq(0)] * (N + 1)
    for k, ck in enumerate(coeffs_u):
        for j in range(k + 1):
            a[j] += ck * comb(k, j) * (-h) ** (k - j) * (2 * h) ** j
    beta = []
    for k in range(N + 1):
        s = fmpq(0)
        for j in range(k + 1):
            s += fmpq(comb(k, j), comb(N, j)) * a[j]
        beta.append(s)
    return min(beta)


def GM80(r):
    M = 80
    a = arb(1) / 4 + M
    b = arb(1) / 4
    t = (r / 2) ** 2
    val = arb.log((a + arb(1) / 2) ** 2 + t) / 2
    val -= a / (a ** 2 + t)
    s = arb(0)
    for k in range(M):
        ak = b + k
        s += ak / (ak ** 2 + t)
    val -= s
    val -= arb(1) / ((a + arb(1) / 2) ** 2 + t) / 4
    val -= arb.log(arb.pi())
    A2 = arb(2).sqrt() * arb(2).log()
    val -= A2 * arb.cos(r * arb(2).log())
    val -= (arb(25) / (2 * arb.pi())).log() - arb(1) / 25 - A2
    return val


# ---------- (2)+(3) per piece ----------
h = fmpq(1, 8)
worst_bern = None
worst_step3 = None
err_max = fmpq(0)
bad2, bad3 = [], []
for (i, c, p_u, _eps) in gen:
    _c2, lb, err = bm.piece_lower(i)
    err_max = max(err_max, err)
    p_r = cd.r_coeffs(c, p_u, fmpq(0))
    d = [lb[n] - p_u[n] for n in range(len(lb))]
    d[0] -= err
    bmin = bernstein_min(d, h)
    if worst_bern is None or bmin < worst_bern[0]:
        worst_bern = (bmin, i)
    if bmin < 0:
        bad2.append(i)
    # step (3): 41 samples per piece
    for j in range(41):
        u = -h + fmpq(j, 20) * h  # -1/8 .. 1/8
        r = c + u
        lbval = sum(lb[n] * u ** n for n in range(len(lb))) - err
        lb_arb = arb(lbval.numerator) / lbval.denominator
        gap = GM80(arb(r.numerator) / r.denominator) - lb_arb
        if worst_step3 is None or float(gap.mid()) < worst_step3[0]:
            worst_step3 = (float(gap.mid()), i, float(r), float(gap.lower()))
        if float(gap.lower()) < 0:
            bad3.append((i, float(r), float(gap.lower())))

print(f"(2) Bernstein certificate p_i <= lb_i - err_i on every piece: "
      f"min Bernstein coeff = {float(worst_bern[0]):.4e} (piece {worst_bern[1]}); "
      f"failures: {bad2 if bad2 else 'none'}")
print(f"    max err_i = {float(err_max):.4e}")
print(f"(3) lb_i - err_i <= G_80 at 4100 samples: min gap = {worst_step3[0]:.4e} "
      f"(piece {worst_step3[1]}, r = {worst_step3[2]}, lower ball {worst_step3[3]:.4e}); "
      f"violations: {bad3[:5] if bad3 else 'none'}")

verdict = "PASS" if (not mismatch and not bad2 and not bad3) else "FAIL"
print("AUDIT:", verdict)
json.dump({
    "table_match": not mismatch, "mismatch_pieces": mismatch,
    "bernstein_min": [float(worst_bern[0]), worst_bern[1]], "bernstein_failures": bad2,
    "err_max": float(err_max),
    "step3_min_gap": worst_step3, "step3_violations": bad3[:20],
    "verdict": verdict,
}, open("/home/hatch/workspace/gppdiscovery3/research/flint_inputb_generators_audit_2026_11_14_results.json", "w"), indent=1)
print("results written")
