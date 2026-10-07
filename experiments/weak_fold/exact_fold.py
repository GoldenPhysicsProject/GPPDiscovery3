"""Exact finite test of a linear exchange fold; not a Standard Model derivation."""
from fractions import Fraction as Q
import json
from pathlib import Path

def zero(n,m): return [[Q(0) for _ in range(m)] for _ in range(n)]
def eye(n): return [[Q(i==j) for j in range(n)] for i in range(n)]
def add(A,B): return [[a+b for a,b in zip(r,s)] for r,s in zip(A,B)]
def scale(c,A): return [[c*a for a in r] for r in A]
def mul(A,B): return [[sum(a*b for a,b in zip(r,c)) for c in zip(*B)] for r in A]
def transpose(A): return [list(r) for r in zip(*A)]
def comm(A,B): return add(mul(A,B),scale(-1,mul(B,A)))
def block(X,k):
    A=zero(4,4)
    for i in range(2):
        for j in range(2): A[i+2*k][j+2*k]=Q(X[i][j])
    return A

I=eye(4)
D=[[Q(j==(i+2)%4) for j in range(4)] for i in range(4)]
P=scale(Q(1,2),add(I,D))
assert mul(D,D)==I and mul(P,P)==P
E=[[0,1],[0,0]]; F=[[0,0],[1,0]]; H=[[1,0],[0,-1]]
for X in [E,F,H]:
    L,R=block(X,0),block(X,1)
    assert mul(mul(D,L),D)==R
    assert mul(mul(P,L),P)==mul(mul(P,R),P)
    assert mul(mul(P,add(L,scale(-1,R))),P)==zero(4,4)
    assert comm(P,add(L,R))==zero(4,4)
    assert comm(P,L)!=zero(4,4)

# On the folded two-dimensional coordinate space, each separate compression
# is X/2. Only the sum acts as X with the original Lie bracket normalization.
EL,FL,HL=[scale(Q(1,2),X) for X in [E,F,H]]
assert comm(EL,FL)==scale(Q(1,2),HL)
assert comm(E,F)==H

# Reproduce Sol's selected incidence Gram; topology does not select L vs R.
B=zero(8,5)
for col,(i,j) in enumerate([(0,1),(1,2),(2,3),(4,5),(6,7)]):
    B[i][col]=Q(1);B[j][col]=Q(-1)
G=mul(transpose(B),B)
expected=[[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,0,0],[0,0,0,2,0],[0,0,0,0,2]]
assert G==expected
permutation=[0,1,2,3,6,7,4,5]
edge_swap=[0,1,2,4,3]
assert [[B[permutation[i]][j] for j in range(5)] for i in range(8)]==[[B[i][edge_swap[j]] for j in range(5)] for i in range(8)]
result={
 'arithmetic':'exact Fraction arithmetic; no numerical tolerance',
 'incidence_gram':expected,
 'weak_block_exchange_is_unweighted_graph_automorphism':True,
 'folded_left_equals_folded_right':True,
 'relative_generators_compress_to_zero':True,
 'diagonal_generators_preserve_fold':True,
 'left_only_generators_do_not_preserve_fold':True,
 'scope':'Linear exchange on two internal doublets only. Does not identify Lorentz chirality, CPT, or an anti-linear charge-conjugation fold.',
 'conclusion':'Symmetric folding selects the diagonal internal action, not one of the two original independent weak factors.'}

# Anti-linear conjugate fold, represented exactly on real coordinates.
# v=(Re v1, Im v1, Re v2, Im v2); C is complex conjugation.
C=[[Q(i==j)*(-1 if i%2 else 1) for j in range(4)] for i in range(4)]
J=zero(4,4)
for i in [0,2]: J[i][i+1]=-Q(1); J[i+1][i]=Q(1)
def diag2(A,B):
    n=len(A);Z=zero(2*n,2*n)
    for i in range(n):
        for j in range(n): Z[i][j]=A[i][j];Z[n+i][n+j]=B[i][j]
    return Z
Dc=zero(8,8)
for i in range(4):
    for j in range(4): Dc[i][j+4]=C[i][j];Dc[i+4][j]=C[i][j]
emb=eye(4)+C
Jphys=diag2(J,scale(-1,J))
assert mul(Dc,emb)==emb
assert mul(Jphys,emb)==mul(emb,J)
assert mul(Jphys,Jphys)==scale(-1,eye(8))
assert comm(Dc,Jphys)==zero(8,8)
# For any complex-linear generator A, its conjugate block is C A C.
# Verify the U(1) generator and a nontrivial real SU(2) generator.
Aweak=[[Q(0),Q(0),Q(1),Q(0)],[Q(0),Q(0),Q(0),Q(1)],
       [Q(-1),Q(0),Q(0),Q(0)],[Q(0),Q(-1),Q(0),Q(0)]]
for A in [J,Aweak]:
    Ad=diag2(A,mul(mul(C,A),C))
    assert mul(Ad,emb)==mul(emb,A)
    assert comm(Ad,Dc)==zero(8,8)
result['conjugate_fold']={
 'fixed_graph':'(v, conjugate(v))',
 'physical_complex_structure':'diag(i,-i)',
 'preserves_original_charged_complex_doublet':True,
 'adds_independent_weak_group':False,
 'limitation':'Lorentz Weyl representation and gauge action must still be supplied; their chirality is not derived.'}
out=Path(__file__).with_name('results.json')
out.write_text(json.dumps(result,indent=2)+'\n')
print(out.read_text())
