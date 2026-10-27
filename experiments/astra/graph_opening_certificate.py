"""Exact graph-opening certificates. No amplitude or duality weights computed."""
from fractions import Fraction
import json

def determinant(a):
    a = [[Fraction(x) for x in row] for row in a]
    result = Fraction(1)
    for k in range(len(a)):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            result = -result
        p = a[k][k]
        result *= p
        for i in range(k+1, len(a)):
            q = a[i][k]/p
            for j in range(k+1, len(a)):
                a[i][j] -= q*a[k][j]
            a[i][k] = 0
    assert result.denominator == 1
    return int(result)

def ladder(L):
    assert L >= 1
    v = 2*(L+1)
    edges = [(2*j, 2*j+1) for j in range(L+1)]
    edges += [(2*j+r, 2*j+2+r) for j in range(L) for r in (0,1)]
    tree = [e for e in edges if e[0]%2 == 0]
    cut = [e for e in edges if e not in tree]
    reached = {0}
    while True:
        new = reached | {b for a,b in tree if a in reached} | {a for a,b in tree if b in reached}
        if new == reached:
            break
        reached = new
    assert len(reached) == v and len(tree) == v-1
    assert len(edges)-v+1 == L == len(cut)
    matrix = [[0]*v for _ in range(v)]
    for a,b in edges:
        matrix[a][a] += 1; matrix[b][b] += 1
        matrix[a][b] -= 1; matrix[b][a] -= 1
    count = determinant([row[:-1] for row in matrix[:-1]])
    return dict(loops=L, vertices=v, internal_edges=len(edges),
                opened_legs=4+2*L, spanning_trees=count,
                chosen_cut_edges=cut)

if __name__ == "__main__":
    results = [ladder(L) for L in range(1,11)]
    counts = [1]+[r["spanning_trees"] for r in results]
    assert counts[1] == 4
    assert all(counts[i] == 4*counts[i-1]-counts[i-2] for i in range(2,len(counts)))
    print(json.dumps(results, indent=2))
