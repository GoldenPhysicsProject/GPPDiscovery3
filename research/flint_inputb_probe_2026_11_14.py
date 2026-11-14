"""
flint_inputb_probe_2026_11_14.py -- cheapest falsifier for input (b) of the
first prime-visible window certificate (GPPVerify PR #597, head c930c5e).

Input (b), as consumed by GppWindowL036.weilForm_autocorr_nonneg_l036:
  (b1) forall r in [-25,25], mfun r <= symbolTwo r - betaTwo 25
  (b2) forall r in [-25,25], |mfun r| <= 8
where (definitions read firsthand from the branch sources):
  symbolTwo r  = Re(digamma(1/4 + i r/2)) - log pi - sqrt(2) log 2 * cos(r log 2)
  betaTwo 25   = log(25/(2 pi)) - 1/25 - sqrt(2) log 2
  mfun r       = m(|r|); on (i/4,(i+1)/4], i<100, m is the degree-12 polynomial
                 pieceTable[i] (ascending rational coefficients,
                 GppWindowMinorantL036, extracted verbatim below).

Method: Arb (python-flint) at prec 200. Piece polynomials are evaluated in
exact rational arithmetic (fmpq) and only then converted to balls, so the
minorant side carries NO rounding error beyond the final conversion.
The digamma side is an acb ball. A grid point is a CERTIFIED violation iff
the upper ball of (mfun) exceeds the lower ball of (symbolTwo - betaTwo).
Grid: every piece sampled at 9 interior points plus both endpoints (endpoint
values taken from BOTH adjacent pieces: the defined value at a breakpoint
comes from the left piece, but a jump up on the right piece is where a
polynomial minorant would fail in practice, so both are reported), plus a
uniform 4001-point grid on [0,25]. mfun is even, so [0,25] covers (b).

Falsifier stated up front: if any sampled point certifiably violates (b1) or
(b2), input (b) as stated is KILLED and the counterexample is the deliverable.
If no violation: report the minimum margin and its location (NUMERICAL,
sampling-based; not a proof).

Also reports the decomposition margins: symbolTwo-betaTwo-G_80 (should be
>= 0 by the banked GM_le, ~0.28/80^2 scale) and G_80 - mfun (the generator's
own slack), to locate where the total margin is spent.

Flint, 2026-11-14. Run: ~/workspace/venvs/arb-venv/bin/python (python-flint).
"""

import json
import re
import sys

import flint
from flint import acb, arb, ctx, fmpq

ctx.prec = 200

LEAN_FILE = "/tmp/WindowMinorantL036.lean"

# ---------- extract pieceTable verbatim from the Lean source ----------
src = open(LEAN_FILE).read()
m = re.search(r"def pieceTable : List \(List ℚ\) := \[(.*?)\n\]", src, re.S)
assert m, "pieceTable not found"
fracs = re.findall(r"\((-?\d+) / (\d+) : ℚ\)", m.group(1))
assert len(fracs) == 100 * 13, f"expected 1300 coefficients, got {len(fracs)}"
pieces = []
for i in range(100):
    coeffs = [fmpq(int(n), int(d)) for (n, d) in fracs[13 * i: 13 * (i + 1)]]
    pieces.append(coeffs)
print(f"pieceTable extracted: 100 pieces x 13 coefficients (exact rationals)")


def mfun_exact(r_q):
    """mfun at a nonnegative rational point, exact fmpq.
    Piece convention: piece i on (i/4, (i+1)/4]; r=0 -> piece 0."""
    if r_q == 0:
        i = 0
    else:
        i = int((r_q * 4 - 1).__ceil__())  # smallest i with r <= (i+1)/4
        i = max(0, min(99, i))
    acc = fmpq(0)
    for c in reversed(pieces[i]):
        acc = acc * r_q + c
    return acc, i


def GM80(r):
    """G_M at M=80, arb ball, exactly as in GppWindowGM.GM."""
    M = 80
    a = arb(1) / 4 + M
    b = arb(1) / 4
    t = (r / 2) ** 2
    val = arb.log(a + arb(1) / 2) if False else None
    # 1/2 log((M+3/4)^2 + (r/2)^2)
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
    val -= (arb(25) / (2 * arb.pi())).log() - arb(1) / 25 - A2  # betaTwo 25
    return val


def target(r):
    """symbolTwo(r) - betaTwo(25), acb/arb balls."""
    z = acb(arb(1) / 4, r / 2)
    re_psi = z.digamma().real
    A2 = arb(2).sqrt() * arb(2).log()
    beta = (arb(25) / (2 * arb.pi())).log() - arb(1) / 25 - A2
    return re_psi - arb.log(arb.pi()) - A2 * arb.cos(r * arb(2).log()) - beta


# ---------- grid ----------
points = set()
for i in range(100):
    for j in range(10):  # endpoints + 8 interior points per piece
        points.add(fmpq(i * 10 + j, 40))  # (i + j/10)/4 in units: r = (10i+j)/40
for n in range(4001):
    points.add(fmpq(n, 160))
points = sorted(p for p in points if p <= 25)

worst_b1 = None       # (margin_mid, r, margin_lower, margin_upper)
worst_decomp_gm = None
worst_decomp_tot = None
max_abs_m = arb(0)
violations = []

for q in points:
    r = arb(q.numerator) / q.denominator
    mv, piece = mfun_exact(q)
    m_ball = arb(mv.numerator) / mv.denominator
    tgt = target(r)
    margin = tgt - m_ball
    lo, hi = float(margin.lower()), float(margin.upper())
    mid = float(margin.mid())
    if worst_b1 is None or mid < worst_b1[0]:
        worst_b1 = (mid, float(q), lo, hi)
    if hi < 0:
        violations.append((float(q), lo, hi))
    if abs(m_ball) > max_abs_m:
        max_abs_m = abs(m_ball)
    gm = GM80(r)
    d_gm = gm - m_ball
    d_tot = tgt - gm
    if worst_decomp_gm is None or float(d_gm.mid()) < worst_decomp_gm[0]:
        worst_decomp_gm = (float(d_gm.mid()), float(q), float(d_gm.lower()))
    if worst_decomp_tot is None or float(d_tot.mid()) < worst_decomp_tot[0]:
        worst_decomp_tot = (float(d_tot.mid()), float(q), float(d_tot.lower()))

print(f"grid points evaluated: {len(points)}")
print(f"(b1) min margin symbolTwo-betaTwo-mfun = {worst_b1[0]:.6e} at r = {worst_b1[1]}")
print(f"     ball enclosure there: [{worst_b1[2]:.6e}, {worst_b1[3]:.6e}]")
print(f"(b2) max |mfun| = {float(max_abs_m.mid()):.6f}  (bound is 8)")
print(f"decomp: min (G_80 - mfun) = {worst_decomp_gm[0]:.6e} at r = {worst_decomp_gm[1]} (lower ball {worst_decomp_gm[2]:.6e})")
print(f"decomp: min (target - G_80) = {worst_decomp_tot[0]:.6e} at r = {worst_decomp_tot[1]} (lower ball {worst_decomp_tot[2]:.6e})")
if violations:
    print(f"*** CERTIFIED VIOLATIONS of (b1) at {len(violations)} sampled points: {violations[:10]}")
    verdict = "KILLED"
else:
    print("no certified violation on the grid: input (b) SURVIVES the probe (NUMERICAL)")
    verdict = "SURVIVES"

json.dump({
    "verdict": verdict,
    "points": len(points),
    "min_margin_b1": worst_b1,
    "max_abs_mfun": float(max_abs_m.mid()),
    "min_G80_minus_mfun": worst_decomp_gm,
    "min_target_minus_G80": worst_decomp_tot,
    "violations": violations[:20],
}, open("/home/hatch/workspace/gppdiscovery3/research/flint_inputb_probe_2026_11_14_results.json", "w"), indent=1)
print("results written to flint_inputb_probe_2026_11_14_results.json")
