import sympy as sp

E, k, mu, s = sp.symbols('E k mu s', positive=True, real=True)

# Explicit CM massive spinors for
# p=(E,0,0,k), q=(E,0,0,-k), with E^2-k^2=mu^2.
# Up to harmless little-group/crossing signs, the angle/square matrices can be
# represented by the following 2x2 matrices.
a2 = E + k
b2 = E - k
A = sp.Matrix([[0, a2], [-b2, 0]])
B = sp.Matrix([[0, b2], [-a2, 0]])
C = A * B

trC = sp.factor(sp.trace(C))
detC = sp.factor(C.det())
sym2 = sp.factor(trC**2 - detC)

assert sp.simplify(trC**2 - 4*(E**2 + k**2)**2) == 0
assert sp.simplify(detC - (E**2-k**2)**2) == 0

# Impose mu^2=E^2-k^2 and s=4E^2.
sym2_mu = sp.expand(sym2.subs(k**2, E**2-mu**2))
expected_E = sp.expand((4*E**2 - 2*mu**2)**2 - mu**4)
assert sp.simplify(sym2_mu - expected_E) == 0
expected_s = sp.expand((s - 2*mu**2)**2 - mu**4)
assert sp.simpl(expected_s - (s**2 - 4*s*mu**2 + 3*mu**4)) == 0

# Threshold k=0, E=mu, s=4 mu^2: three vector polarizations.
threshold = sp.simplify(expected_s.subs(s, 4*mu**2))
assert threshold == 3*mu**4

print('tr(C) =', trC)
print('det(C) =', detC)
print('Tr_Sym2(C) =', sym2)
print('Invariant form = (s - 2 mu^2)^2 - mu^4')
print('Threshold =', threshold)
