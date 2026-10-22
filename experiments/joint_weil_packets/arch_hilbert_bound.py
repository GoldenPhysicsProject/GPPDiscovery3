"""Finite numerical checks of the exact archimedean Hankel bound; not a Lean certificate."""
import numpy as np
from scipy.special import polygamma
from scipy.linalg import eigvalsh
import json
rows=[]
for N in [8,32,128,512]:
 j=np.arange(N)
 G=1/(2*j[:,None]+2*j[None,:]+1)
 lam=float(eigvalsh(G,subset_by_index=[N-1,N-1])[0])
 assert lam<np.pi/2
 rows.append({'rank':N,'largest_gram_eigenvalue':lam,'upper_bound':np.pi/2})
# Exact-series quadratic form on exp(-a*u): (1/4)*trigamma((a+1/2)/2).
probes=[]
for a in [.1,.5,1,5,100]:
 q=.25*polygamma(1,(a+.5)/2)
 ratio=float(2*a*q)
 assert ratio<np.pi/2
 probes.append({'a':a,'quadratic_form_over_L2norm_squared':ratio})
print(json.dumps({'scope':__doc__,'gram_checks':rows,'exponential_probes':probes},indent=2))
