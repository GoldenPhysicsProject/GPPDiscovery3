"""Finite CAR realization of paired orientations; no zeta-zero identification."""
import numpy as np
I=np.eye(2); Z=np.diag([1.,-1.]); a=np.array([[0.,1.],[0.,0.]])
c=[np.kron(a,I),np.kron(Z,a)]; eye=np.eye(4)
n=[x.T@x for x in c]; H=2*(n[0]+n[1]); K=2*(n[0]-n[1]); Q=n[0]-n[1]
# Fermionic swap: double occupation acquires a minus sign.
S=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,-1.]])
assert np.array_equal(S@S,eye)
for i in range(2):
 for j in range(2):
  assert np.array_equal(c[i]@c[j]+c[j]@c[i],np.zeros((4,4)))
  assert np.array_equal(c[i]@c[j].T+c[j].T@c[i],eye if i==j else np.zeros((4,4)))
 assert np.array_equal(S@c[i]@S,c[1-i])
assert np.array_equal(S@H@S,H)
assert np.array_equal(S@K@S,-K)
assert np.array_equal(S@Q@S,-Q)
assert np.linalg.eigvalsh(H).min()==0
b=(c[0]+c[1])/np.sqrt(2)
assert np.allclose(b@b.T+b.T@b,eye)
assert np.allclose(S@b@S,b)
d=(c[0]-c[1])/np.sqrt(2)
assert np.allclose(K@b-b@K,-2*d)
assert np.allclose(K@d-d@K,-2*b)
assert np.allclose(H@b-b@H,-2*b)
# Two signed branches with t_minus=-tau have the same local positive-energy phase.
for tau in [-2.,.3,1.7]:
 assert abs(np.exp(-1j*2*tau)-np.exp(-1j*(-2)*(-tau)))<1e-14
print('PASS: two-mode CAR; graded swap; H positive/invariant; K,Q odd; normalized shared CAR mode; local phase identity.')
print('H eigenvalues:',np.linalg.eigvalsh(H).tolist())
