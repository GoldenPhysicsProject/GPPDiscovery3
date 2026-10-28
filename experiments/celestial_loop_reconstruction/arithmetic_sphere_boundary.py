"""Finite interacting CAR boundary on S2 and Poisson extension. Not a celestial CFT."""
import math,itertools,json
primes=(2,3)
modes=[(p,l,a) for p in primes for l,aa in [(0,range(1)),(1,range(3))] for a in aa]
states=list(itertools.product((0,1),repeat=len(modes)))
nodes=[tuple(s if k==j else 0 for k in range(3)) for j in range(3) for s in (-1,1)]
nodes=[tuple(n) for n in nodes]
def basis(n):return [1/math.sqrt(4*math.pi)]+[math.sqrt(3/(4*math.pi))*x for x in n]
# Exact quadrature through degree 2 on the octahedron.
for a in range(4):
 for b in range(4):
  val=sum(basis(n)[a]*basis(n)[b] for n in nodes)*4*math.pi/6
  assert abs(val-(a==b))<1e-14
# CAR with ordering signs, explicitly tested on every occupation state.
def op(state,i,create):
 if state[i]==int(create):return None
 t=list(state);t[i]=int(create)
 return (-1)**sum(state[:i]),tuple(t)
def compose(state,i,ci,j,cj):
 first=op(state,j,cj)
 if first is None:return {}
 second=op(first[1],i,ci)
 return {} if second is None else {second[1]:first[0]*second[0]}
for state in states:
 for i in range(8):
  for j in range(8):
   for cj in [False,True]:
    d=compose(state,i,False,j,cj)
    for s,v in compose(state,j,cj,i,False).items():d[s]=d.get(s,0)+v
    d={s:v for s,v in d.items() if v}
    assert d==({state:1} if cj and i==j else {})
def calculate(g,beta=1.3):
 def energy(s):
  # Angular kinetic coefficient chosen as 0.2; not derived from arithmetic.
  free=sum(n*(math.log(p)+.2*l*(l+1)) for n,(p,l,a) in zip(s,modes))
  return free+g*sum(s[:4])*sum(s[4:])
 es=[energy(s) for s in states];ws=[math.exp(-beta*e) for e in es];Z=sum(ws)
 def mean(f):return sum(w*f(s) for s,w in zip(states,ws))/Z
 cov=mean(lambda s:sum(s[:4])*sum(s[4:]))-mean(lambda s:sum(s[:4]))*mean(lambda s:sum(s[4:]))
 # Rotational invariance makes occupations within each l=1 multiplet equal.
 for start in [0,4]:
  occ=[mean(lambda s,i=i:s[i]) for i in range(start+1,start+4)]
  assert max(occ)-min(occ)<1e-14
 if g==0:
  expected=math.prod(1+math.exp(-beta*(math.log(p)+.2*l*(l+1))) for p,l,a in modes)
  assert abs(Z-expected)<1e-12 and abs(cov)<1e-13
 else:assert abs(cov)>1e-3
 return dict(g=g,partition=Z,minimum_energy=min(es),connected_prime_number_four_point=cov)
# Harmonic solid-ball extension: l=0 is constant; l=1 is linear in Cartesian x.
# Boundary reconstruction of coefficients via the six-node quadrature.
a=[.3,-.2,.7,1.1]
f=[sum(x*y for x,y in zip(a,basis(n))) for n in nodes]
back=[4*math.pi/6*sum(v*basis(n)[j] for v,n in zip(f,nodes)) for j in range(4)]
assert max(abs(x-y) for x,y in zip(a,back))<1e-14
result={'CAR':'all 256 states and all mode pairs passed','angular_gram_and_inverse':'passed','models':[calculate(0),calculate(.4)],'scope':'Chosen finite interacting S2 model; harmonic ball extension, not Lorentzian holography or derived gravity'}
print(json.dumps(result,indent=2))
open('arithmetic_sphere_boundary_results.json','w').write(json.dumps(result,indent=2)+'\n')
