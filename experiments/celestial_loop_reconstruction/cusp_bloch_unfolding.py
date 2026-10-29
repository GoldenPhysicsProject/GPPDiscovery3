import numpy as np,json
from scipy.integrate import simpson
x=np.linspace(0,1,257)
th=2*np.pi*np.arange(128)/128
n=np.arange(-8,9)
f=lambda q:np.pi**(-.25)*np.exp(-q*q/2)
values=f(x[:,None]+n[None,:])
B=values@np.exp(-1j*n[:,None]*th[None,:])
fiber_norm=np.mean(abs(B)**2,axis=1)
norm1=float(simpson(fiber_norm,x=x))
errs=[]
for j in [-3,-1,0,2,3]:
 reconstructed=np.mean(B*np.exp(1j*j*th)[None,:],axis=1)
 errs.append(float(np.max(abs(reconstructed-f(x+j)))))
Bshift=f(x[:,None]+1+n[None,:])@np.exp(-1j*n[:,None]*th[None,:])
quasi=float(np.max(abs(Bshift-np.exp(1j*th)[None,:]*B)))
assert abs(norm1-1)<1e-12 and max(errs)<1e-12 and quasi<1e-10
# Full two-dimensional Gaussian is separable, so its Bloch norm is norm1 squared.
result={'one_dimensional_bloch_norm':norm1,'two_dimensional_bloch_norm':norm1**2,'inverse_reconstruction_max_error':max(errs),'quasiperiodic_boundary_max_error':quasi,'repeat_cell_norm_squared':{str(N):N for N in [1,4,16,64]},'normalized_repeat_mass_in_one_cell':{str(N):1/N for N in [1,4,16,64]},'twisted_cusp_transverse_gap_examples':{str(t):t*t for t in [0,.1,.5,1.]},'scope':'Finite Gaussian Bloch quadrature, not full nonabelian Bianchi spectral decomposition.'}
print(json.dumps(result,indent=2))
