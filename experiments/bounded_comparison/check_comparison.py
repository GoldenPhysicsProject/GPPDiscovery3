"""Exact rational tests of the finite bounded-comparison identities. Not RH."""
import sympy as s
M=s.Matrix([[4,1,1],[1,3,0],[1,0,2]])
E=s.Matrix([[1,0],[0,1],[0,0]])
m=s.Rational(1,2)
I=s.eye(2)
A=E.T*M*E
B=E.T*M*M*E
G=E.T*M.inv()*E
C=B/m-A
D=A/m-I
Xs=C.inv()*D

def W(X): return I/m+X.T*C*X-D*X-X.T*D

def zero(X): assert X==s.zeros(*X.shape), X

assert all((M-m*s.eye(3))[:k,:k].det()>0 for k in range(1,4))
assert all(C[:k,:k].det()>0 for k in range(1,3))
for X in [s.zeros(2),s.eye(2),s.Matrix([[1,2],[-3,4]]),Xs]:
    R=E-M*E*X
    zero(W(X)-G-R.T*(s.eye(3)/m-M.inv())*R)
    zero(W(X)-W(Xs)-(X-Xs).T*C*(X-Xs))
zero(B-A*A-E.T*M*(s.eye(3)-E*E.T)*M*E)
assert W(Xs)-G!=s.zeros(2)
assert A*B!=B*A # Noncommutative test, not a scalar surrogate.
print('PASS: 8 exact comparison identities, variance identity, strict-positive noncommuting fixture.')
print('Optimal residual is nonzero:', W(Xs)-G)
print('Scope: rational finite algebra only; no analytic quadrature or RH claim.')
