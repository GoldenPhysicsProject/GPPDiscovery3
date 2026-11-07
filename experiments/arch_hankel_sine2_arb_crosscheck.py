#!/usr/bin/env python3
"""Independent Arb check of the EXACT rational Wave234 cross-Arch intervals.

Optional dependency: python-flint (this is Grok's interval-arithmetic style).
Example:
  python3 arch_hankel_sine2_k6_exact_intervals.py | \
      python3 arch_hankel_sine2_arb_crosscheck.py

The exact rational interval producer remains the authoritative proof;
Arb independently checks that its 16-term head is inside all enclosures.
No attempt to prove the tail using this floating ball comparison.
"""

import json
import sys

from flint import arb, ctx, fmpq

ctx.prec = 256


def main() -> None:
    data = json.load(sys.stdin)
    scale = int(data["scale"])
    pi2 = arb.pi() ** 2
    positions = [arb(fmpq(7+14*j, 20)) for j in range(6)]
    laps = []
    for n in range(16):
        lam = arb(fmpq(4*n+1, 2))
        s = 2*pi2*(1-(-lam).exp()) / (lam*(lam*lam+4*pi2))
        laps.append([(-lam*a).exp()*s for a in positions])

    checked = 0
    for j in range(6):
        for k in range(6):
            head = -sum((v[j]*v[k] for v in laps), arb(0))
            entry = data["matrix"][j][k]
            lo = arb(fmpq(int(entry["lower_scaled"]), scale))
            hi = arb(fmpq(int(entry["upper_scaled"]), scale))
            assert head > lo, (j, k, "below lower", str(head), str(lo))
            assert head < hi, (j, k, "above upper", str(head), str(hi))
            checked += 1
    print(f"Arb 256-bit independent cross-check: {checked}/36 head entries "
          "enclosed by rational bounds; PASS")


if __name__ == "__main__":
    main()
