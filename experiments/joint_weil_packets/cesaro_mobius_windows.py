"""Exact rational audit of Cesaro/window energy identities, including all boundaries."""
from fractions import Fraction as Q
import json

def audit(a,H):
    r=len(a);M=[0]
    for v in a:M.append(M[-1]+v)
    pair=sum((a[m-1]*a[n-1]*max(Q(0),1-Q(abs(m-n),H)) for m in range(1,r+1) for n in range(1,r+1)),Q(0))
    windows=sum((sum(a[n-1] for n in range(1,r+1) if j<=n<j+H)**2 for j in range(2-H,r+1)),Q(0))/H
    if H<=r:
        boundary=(sum(M[k]**2 for k in range(1,H))+sum((M[r]-M[r-k])**2 for k in range(1,H)))/Q(H)
        interior=sum((M[j+H-1]-M[j-1])**2 for j in range(1,r-H+2))/Q(H)
    else:
        boundary=sum(M[k]**2+(M[r]-M[r-k])**2 for k in range(1,r))/Q(H)
        interior=Q(H-r+1,H)*M[r]**2
    assert pair==windows==boundary+interior
    if H>=r:assert windows>=Q(H-r+1,H)*M[r]**2
    return windows,boundary,interior

mu=[0,1,-1,-1,0,-1,1,-1,0,0,1,-1,0,-1,1,1]
checks=[]
for N in [3,8,16]:
    a=mu[1:N];M=[0]
    for x in a:M.append(M[-1]+x)
    E=sum((Q(M[r]**2,r*(r+1)) for r in range(1,N)),Q(0))
    D=sum((a[n-1]**2*(Q(1,n)-Q(1,N)) for n in range(1,N)),Q(0))
    for H in [1,2,3,N,2*N]:
        total=sum((audit(a[:r],H)[0]/(r*(r+1)) for r in range(1,N)),Q(0))
        direct=sum((a[m-1]*a[n-1]*(Q(1,max(m,n))-Q(1,N))*max(Q(0),1-Q(abs(m-n),H)) for m in range(1,N) for n in range(1,N)),Q(0))
        assert total==direct and 0<=total<=H*D
        if H>=N-1:assert total>=Q(H-N+2,H)*E
        checks.append({'N':N,'H':H,'Cesaro_energy':str(total),'original_energy':str(E),'exact':True})
# Constant-sign control: same positive identities hold, but energy is of order N.
for N in [8,16]:
    a=[1]*(N-1)
    for H in [1,3,N,2*N]:audit(a,H)
print(json.dumps({'scope':__doc__,'checks':checks,'positive_control_exact':True},indent=2))
