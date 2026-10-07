"""Exact integer checks of observable descent; no physical CPT claim."""
import numpy as np
D = np.fliplr(np.eye(4, dtype=int))
Rq = np.eye(4, dtype=int)[[2, 3, 0, 1]]
Rt = np.eye(4, dtype=int)[[1, 0, 3, 2]]
assert np.array_equal(Rq @ D, Rt)
# Symmetrized matrix units span the invariant observable algebra.
for i in range(4):
    for j in range(4):
        E = np.zeros((4, 4), dtype=int)
        E[i, j] = 1
        O = E + D @ E @ D
        assert np.array_equal(Rq @ O @ Rq, Rt @ O @ Rt)
v = np.array([1, 0, 0, -1])
assert np.array_equal(D @ v, -v)
assert not np.array_equal(Rq @ v, Rt @ v)
# Twice the twirled density matrix avoids fractions.
rho2 = np.diag([1, 0, 0, 1])
assert np.array_equal(D @ rho2 @ D, rho2)
assert not np.array_equal(D @ rho2, rho2)
print("PASS: observable actions agree; vector equality and fixed support do not follow.")
