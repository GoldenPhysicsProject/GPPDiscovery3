import math,json
import numpy as np
from scipy.integrate import quad
rows=[]
for n in range(5):
 h=2.**(-n)
 # x times the normalized one-dimensional Gaussian is the dipole.
 def err(x):
  g=lambda y: math.pi**(-.25)*math.exp(-y*y/2)
  return ((g(x-h)-g(x+h))/(2*h)-x*g(x))**2
 measured=quad(err,-12,12,epsabs=1e-13)[0]
 expected=.5-math.expm1(-h*h)/(2*h*h)-math.exp(-h*h/4)
 assert abs(measured-expected)<2e-13
 rows.append({'h':h,'dipole_error_squared':measured,'analytic':expected})
# Matrix group checks: det-one arithmetic generators; D acts by z->4z.
D=np.diag([2.,.5]).astype(complex)
W=np.array([[0,-1],[1,0]],complex)
N=lambda a:np.array([[1,a],[0,1]],complex)
res=[]
for n in range(7):
 for a in (1,1j):
  lhs=np.linalg.matrix_power(D,-n)@N(a)@np.linalg.matrix_power(D,n)
  res.append(float(np.max(abs(lhs-N(a/4**n)))))
assert max(res)==0
centers=np.array([0,1,1j,1+1j,.5+.25j,-.5j])
G=np.exp(-abs(centers[:,None]-centers[None,:])**2/4)
eig=np.linalg.eigvalsh(G)
assert eig.min()>0
result={'dipole_approximation':rows,'arithmetic_generator_max_error':max(res),'finite_gaussian_gram_min_eigenvalue':float(eig.min()),'scope':'Numerical examples only; density and Lorentz closure proved analytically in accompanying note.'}
print(json.dumps(result,indent=2))
