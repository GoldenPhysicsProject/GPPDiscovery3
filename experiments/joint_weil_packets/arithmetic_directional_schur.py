"""Prime-side parity matrices and fixed-border variational prices; numerical, not certified."""
import contextlib,io,json,math
with contextlib.redirect_stdout(io.StringIO()):
 exec(open('rh_audit/detector_full_line.py').read().split('zeros=[]')[0])
m.mp.dps=45
K=7
N=math.ceil(math.exp(1.7+1.4*(K-1)+1))
sieve=bytearray(b'\x01')*(N+1);sieve[0:2]=b'\x00\x00'
for p in range(2,math.isqrt(N)+1):
 if sieve[p]:sieve[p*p:N+1:p]=b'\x00'*(((N-p*p)//p)+1)
d.N=N;d.lam=[m.mpf(0)]*(N+1)
for p in range(2,N+1):
 if sieve[p]:
  j=p
  while j<=N:d.lam[j]=m.log(p);j*=p
cache={}
def cc(u):
 key=str(u)
 if key not in cache:cache[key]=C(u)
 return cache[key]
a=[m.mpf('.35')+m.mpf('.7')*j for j in range(K)]
Q=m.matrix([[cc(abs(a[i]-a[j])) for j in range(K)] for i in range(K)])
B=m.matrix([[cc(a[i]+a[j]+1) for j in range(K)] for i in range(K)])
rows=[]
for sign in [-1,1]:
 P=Q+sign*B
 eigmins=[]
 for k in range(2,K+1):
  A=P[:k,:k]
  pivot=A[k-1,k-1]-(A[:k-1,k-1].T*m.lu_solve(A[:k-1,:k-1],A[:k-1,k-1]))[0]
  eigmins.append(dict(K=k,min_eigenvalue=str(m.eigsy(A,eigvals_only=True)[0]),schur_pivot=str(pivot)))
 b=P[:K-1,K-1];d0=P[K-1,K-1]
 prices=[]
 for n in range(1,K):
  A=P[:n,:n];bn=b[:n,:];x=m.lu_solve(A,bn);q=(bn.T*x)[0]
  # Independently price residual in the entire known 6-dimensional prefix.
  xf=m.matrix(K-1,1)
  for j in range(n):xf[j]=x[j]
  Af=P[:K-1,:K-1];r=b-Af*xf
  residual=(r.T*m.lu_solve(Af,r))[0]
  full=(b.T*m.lu_solve(Af,b))[0]
  assert abs(q+residual-full)<m.mpf('1e-35')
  prices.append(dict(n=n,q=str(q),border_minus_q=str(d0-q),known_residual_energy=str(residual)))
 rows.append(dict(parity=sign,margins=eigmins,border_diagonal=str(d0),fixed_border_prices=prices))
print(json.dumps(dict(scope="Prime-only reconstruction of fixed-grid parity forms; no zero data; floating diagnostics, not interval certificates",prime_limit=N,rows=rows),indent=2))
