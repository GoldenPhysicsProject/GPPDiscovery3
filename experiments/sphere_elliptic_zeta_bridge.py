#!/usr/bin/env python3
"""S^2 -> P^1 -> elliptic Hasse-Weil bridge; exact finite arithmetic checks.

No RH/BSD assumption or numerical zero input. BAD is an abstract reciprocal
matrix, explicitly NOT geometric Frobenius of an elliptic curve.
"""
from fractions import Fraction as Q
import json

def point_count_prime(p):
    """E: y^2=x^3-x over F_p, p odd, with point at infinity."""
    return 1 + sum((y*y-x*x*x+x) % p == 0
                   for x in range(p) for y in range(p))

def least_nonsquare(p):
    return next(d for d in range(2,p)
                if pow(d,(p-1)//2,p)==p-1)

def point_count_quadratic(p):
    """Brute force over F_p[t]/(t^2-d), d a non-square."""
    d=least_nonsquare(p)
    pairs=[(a,b) for a in range(p) for b in range(p)]
    def mul(x,y):
        a,b=x; c,e=y
        return ((a*c+d*b*e)%p,(a*e+b*c)%p)
    square_counts={}
    for y in pairs:
        y2=mul(y,y)
        square_counts[y2]=square_counts.get(y2,0)+1
    count=1
    for x in pairs:
        x3=mul(mul(x,x),x)
        rhs=((x3[0]-x[0])%p,(x3[1]-x[1])%p)
        count+=square_counts.get(rhs,0)
    return count

def mul(A,B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))
def transpose(A):
    return tuple(tuple(A[j][i] for j in range(2)) for i in range(2))
def scale(k,A):
    return tuple(tuple(k*x for x in row) for row in A)
def det(A):
    return A[0][0]*A[1][1]-A[0][1]*A[1][0]

def mirror_vs_height():
    q=Q(16)
    J=((Q(0),Q(1)),(Q(-1),Q(0)))
    F_bad=((Q(8),Q(0)),(Q(0),Q(2)))
    F_good=((Q(0),Q(-16)),(Q(1),Q(0)))
    H_good=((Q(1),Q(0)),(Q(0),Q(16)))
    assert all(mul(mul(transpose(F),J),F)==scale(q,J)
               for F in (F_bad,F_good))
    assert det(F_bad)==det(F_good)==q
    assert mul(mul(transpose(F_good),H_good),F_good)==scale(q,H_good)
    assert H_good[0][0]>0 and det(H_good)>0
    # BAD-invariant H=[[a,b],[b,c]] requires (64-16)a=(4-16)c=0;
    # thus a=c=0 and det(H)=-b^2, never positive definite.
    for b in range(-4,5):
        H=((Q(0),Q(b)),(Q(b),Q(0)))
        assert mul(mul(transpose(F_bad),H),F_bad)==scale(q,H)
        assert det(H)<=0
    # Reciprocal characteristic polynomial 1-tr(F)T+qT^2:
    assert (Q(1)-10*Q(1,8)+q*Q(1,8)**2)==0
    assert (Q(1)-10*Q(1,2)+q*Q(1,2)**2)==0
    return {"q":16, "both_symplectic_dual":True,
            "bad_eigenvalues":[2,8],"good_eigenvalue_modulus":4,
            "bad_has_positive_invariant_height":False,
            "good_positive_height_diag":[1,16]}

def verify():
    rows=[]
    for p in (3,5,7,11):
        n1=point_count_prime(p)
        a=p+1-n1
        n2=point_count_quadratic(p)
        assert n2==p*p+1-a*a+2*p
        # For each n, log Z_{E x P1}(T) coefficient is
        # N_E(n)*(1+p^n)/n = coefficient of log Z_E(T)+log Z_E(pT).
        for n,ne in ((1,n1),(2,n2)):
            assert Q(ne*(1+p**n),n)==Q(ne,n)+Q(ne*p**n,n)
        rows.append({"p":p, "points_E_Fp":n1, "points_E_Fp2":n2,
                     "a_p":a, "middle_polynomial":[1,-a,p],
                     "points_P1_Fp":p+1,"points_P1_Fp2":p*p+1})
    return {"curve":"y^2=x^3-x","good_primes":rows,
            "projective_line_over_Z":"zeta(P1_Z,s)=zeta(s)*zeta(s-1)",
            "projective_bundle":"zeta(P1_X,s)=zeta(X,s)*zeta(X,s-1)",
            "local_elliptic":"Z(E/Fp,T)=(1-a_p*T+p*T^2)/((1-T)*(1-p*T))",
            "mirror_falsifier":mirror_vs_height(),
            "limits":"finite checks + analytic identities; not RH, BSD or an arithmetic realization of F_bad"}

if __name__=="__main__":
    print(json.dumps(verify(),sort_keys=True,indent=2))
