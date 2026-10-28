"""Rest Dirac spin and CAR-preserving reconstruction checks; hbar=c=1."""
import numpy as np
from scipy.linalg import expm
I=np.eye(2); X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1.,-1.]); P=[X,Y,Z]
H=np.kron(X,I); J=[np.kron(I,p)/2 for p in P]
assert np.allclose(sum(j@j for j in J),.75*np.eye(4))
assert all(np.allclose(H@j,j@H) for j in J)
for j in J:
 assert np.allclose(expm(-2j*np.pi*j),-np.eye(4))
 assert np.allclose(expm(-4j*np.pi*j),np.eye(4))
W=np.vstack([I,I])/np.sqrt(2)
assert np.allclose(W.conj().T@W,I)
assert np.allclose(H@W,W)
for j,p in zip(J,P):assert np.allclose(j@W,W@(p/2))
# Arithmetic orbital reorganization on four binary prime-labelled modes.
# A chosen unitary demonstrates preservation; its physical identification is NOT derived.
V=np.array([[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]])/2
K=np.kron(V,W)
assert np.allclose(K.conj().T@K,np.eye(8))
for p,j in zip(P,J):
 assert np.allclose(np.kron(np.eye(4),j)@K,K@np.kron(np.eye(4),p/2))
# For CAR c_i, {sum_i conjugate(K_ia)c_i, adjoint_b}=(K†K)_ab.
print('PASS: spin Casimir 3/4, mass commutes all rotations, 2pi/4pi, normalized positive-rest embedding, 4-orbital spin/CAR isometry.')
print('Scope: chosen finite isometry, not a derived arithmetic-to-celestial dictionary.')
