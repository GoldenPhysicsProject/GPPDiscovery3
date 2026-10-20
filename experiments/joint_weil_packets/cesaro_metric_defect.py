"""Finite Cesaro metric checks; model generators, not the zeta zero space."""
import json
import mpmath as mp
mp.mp.dps=60

def avg(A,R):
 return mp.matrix([[mp.quad(lambda t:(mp.expm(t*A).T*mp.expm(t*A))[i,j],[-R,0,R])/(2*R) for j in range(2)] for i in range(2)])
models={'rotation':mp.matrix([[0,-1],[1,0]]),'growth_pair':mp.matrix([[mp.mpf(1)/4,0],[0,-mp.mpf(1)/4]]),'jordan':mp.matrix([[0,1],[0,0]])}
rows=[]
for name,A in models.items():
 for R in map(mp.mpf,[1,4,10]):
  G=avg(A,R); D=A.T*G+G*A
  endpoint=(mp.expm(R*A).T*mp.expm(R*A)-mp.expm(-R*A).T*mp.expm(-R*A))/(2*R)
  err=mp.norm(D-endpoint)
  assert err<mp.mpf('1e-50')
  # All these model averages are diagonal.
  invsqrt=mp.diag([1/mp.sqrt(G[i,i]) for i in range(2)])
  eigen=mp.eigsy(invsqrt*D*invsqrt,eigvals_only=True)
  eta=max(abs(x) for x in eigen)
  if name=='rotation':assert eta<mp.mpf('1e-50')
  if name=='growth_pair':assert abs(eta-mp.mpf('0.5'))<mp.mpf('1e-50')
  if name=='jordan':assert abs(eta-1/mp.sqrt(1+R**2/3))<mp.mpf('1e-50')
  rows.append(dict(model=name,R=str(R),G_diagonal=[mp.nstr(G[i,i],20) for i in range(2)],relative_defect=mp.nstr(eta,20),condition_number=mp.nstr(max(G[0,0],G[1,1])/min(G[0,0],G[1,1]),20),endpoint_identity_error=mp.nstr(err,5)))
print(json.dumps(dict(scope='Finite model identities checked numerically; no arithmetic embedding, interval or Lean proof',rows=rows),indent=2))
