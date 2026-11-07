#!/usr/bin/env python3
"""Merge matrix results and fit the decay law."""
import glob, json, math, sys
import numpy as np

rows = []
for fn in sorted(glob.glob("artifacts/**/result.json", recursive=True)):
    rows.append(json.load(open(fn)))
rows.sort(key=lambda r: (r["N"], r["c"]))
with open("results.jsonl", "w") as f:
    for r in rows:
        f.write(json.dumps(r) + "\n")

G1 = 14.1347251417346937
out = []
out.append("# Truncated Weil ground-state scan\n")
out.append("## N-convergence (fixed c=19, the far end)\n")
out.append("| N | dim | lambda_1 | lambda_2 | gap | log10(g) | odd-in-gap | sec |")
out.append("|---|-----|----------|----------|-----|----------|------------|-----|")
for r in [r for r in rows if r["c"] == 19]:
    out.append(f"| {r['N']} | {r['dim']} | `{r['lam1']}` | `{r['lam2']}` | `{r['gap']}` "
               f"| {r['gap_log10']:.4f} | {r['odd_inside_gap']} | {r['seconds']} |")

best = max((r["N"] for r in rows), default=None)
scan = [r for r in rows if r["N"] == best]
out.append(f"\n## c-scan at N={best}\n")
out.append("| c | a=log c | lambda_1 | log10(l1) | gap | log10(gap) | l1/gap |")
out.append("|---|---------|----------|-----------|-----|------------|--------|")
for r in scan:
    out.append(f"| {r['c']} | {r['a']:.5f} | `{r['lam1']}` | {r['lam1_log10']:.4f} "
               f"| `{r['gap']}` | {r['gap_log10']:.4f} | {r['lam1_log10']-r['gap_log10']:+.4f} |")

pos = [r for r in scan if r["sign1"] > 0]
if len(pos) >= 3:
    A = np.array([r["a"] for r in pos])
    Y = np.array([r["lam1_log10"] for r in pos])
    s, i = np.polyfit(A, Y, 1)
    resid = Y - (s * A + i)
    nat = s * math.log(10)
    out.append(f"\n## Fit (positive points only, n={len(pos)})\n")
    out.append(f"- `log10(lambda) = {s:.6f} * a + {i:.6f}`")
    out.append(f"- natural-log slope: **{nat:.4f}**")
    out.append(f"- `-2*gamma_1 = {-2*G1:.4f}`   ratio {nat/(-2*G1):.4f}")
    out.append(f"- `-4*gamma_1 = {-4*G1:.4f}`   ratio {nat/(-4*G1):.4f}")
    out.append(f"- max |residual| in log10: **{np.abs(resid).max():.4f}**")
    out.append("\n### consecutive slopes (curvature check)\n")
    for j in range(len(pos)-1):
        da=pos[j+1]["a"]-pos[j]["a"]
        dy=pos[j+1]["lam1_log10"]-pos[j]["lam1_log10"]
        out.append(f"- c={pos[j]['c']} -> {pos[j+1]['c']}: slope {dy/da:.4f} (log10 per unit a)")
    out.append("\nIf these are not roughly constant there is no single decay "
               "constant and the power-law picture is wrong.\n")
    out.append("\nA small residual means the decay really is a power law in c "
               "and the constant is meaningful. A large one means there is "
               "curvature and the single-constant picture is wrong.")
open("RESULTS.md", "w").write("\n".join(out) + "\n")
print("\n".join(out))
