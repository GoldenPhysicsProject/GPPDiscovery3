"""Finite Gram checks and cutoff growth; infinite divergence proved separately."""
import itertools, json
import mpmath as mp
mp.mp.dps=50
logs=[mp.log(p) for p in [2,3,5]]
def energies(N,L):
 return [sum((k*a for k,a in zip(v,logs)),mp.mpf(0))**2 for v in itertools.product(range(-L,L+1),repeat=N)]
def trace(t,es):return mp.fsum(mp.exp(-mp.mpf(t)*e) for e in es)
es=energies(3,4)
times=list(map(mp.mpf,['.1','.2','.3','.5','.8']))
M=mp.matrix([[trace(a+b,es) for b in times] for a in times])
eigs=mp.eigsy(M,eigvals_only=True)
rows=[]
for L in [4,8,16,32]:
 e=energies(2,L)
 rows.append(dict(N=2,L=L,trace_t1=str(trace(1,e))))
result=dict(scope=__doc__,three_prime_cutoff4=[dict(t=str(t),trace=str(trace(t,es))) for t in times],eigenvalues=[str(x) for x in eigs],matrix_trace=str(sum(M[i,i] for i in range(5))),cutoff_growth=rows,claimed_eigenvalue_sum=str(sum(map(mp.mpf,['.00173195729','.0682183466','1.80086926','28.9145334','569.522346']))),claimed_matrix_trace=str(sum(map(mp.mpf,['191.6554','137.671','113.0665','88.0466','69.8683']))))
print(json.dumps(result,indent=2))
