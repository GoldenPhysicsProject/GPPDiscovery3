"""Numerical shift-band diagnostic for the exact Mobius energy kernel; no interval certificate."""
import numpy as np
from scipy.signal import correlate
import json
NMAX=1000000
mu=np.ones(NMAX,dtype=np.int8);mu[0]=0
prime=np.ones(NMAX,dtype=bool);prime[:2]=False
for p in range(2,NMAX):
    if prime[p]:
        mu[p::p]*=-1;prime[2*p::p]=False;mu[p*p::p*p]=0
rows=[]
for N in [100,10000,1000000]:
    u=mu[1:N].astype(float);ns=np.arange(1,N,dtype=float)
    v=u*(1/ns-1/N)
    cross=2*correlate(v,u,mode='full',method='fft')[len(u):]
    M=np.cumsum(u);E=float(np.sum(M*M/(ns*(ns+1))))
    D=float(np.dot(u,v));C=float(np.sum(cross))
    assert abs(C-(E-D))<1e-9
    max_direct_error=None
    if N==100:
        direct=np.array([2*np.dot(u[:-h],v[h:]) for h in range(1,len(u))])
        max_direct_error=float(np.max(np.abs(cross-direct)))
        assert max_direct_error<1e-12
    bands=[];lo=1
    while lo<N-1:
        hi=min(2*lo,N-1); vals=cross[lo-1:hi-1]
        bands.append({'h_from':lo,'h_below':hi,'signed':float(np.sum(vals)),'absolute':float(np.sum(abs(vals)))})
        lo=hi
    rows.append({'N':N,'signed_cross':C,'absolute_shift_sum':float(np.sum(abs(cross))),'direct_check_error':max_direct_error,'bands':bands})
print(json.dumps({'scope':__doc__,'results':rows},indent=2))
