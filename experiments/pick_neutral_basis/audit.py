"""Decimal Euler--Maclaurin diagnostics for a pole-neutral Pick basis."""
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
import json
import math
from pathlib import Path

GAMMA = D('0.5772156649015328606065120900824024310421593359399235988057672348848677')
PI = D('3.1415926535897932384626433832795028841971693993751058209749445923078164')

def dec(x):
    return D(x.numerator)/D(x.denominator)

def bernoulli(n):
    a = []
    out = []
    for m in range(n+1):
        a.append(F(1,m+1))
        for j in range(m,0,-1):
            a[j-1] = j*(a[j-1]-a[j])
        out.append(a[0])
    return out

B = bernoulli(32)

def zeta_pair(s, n):
    ss, nn = D(s), D(n)
    ln = nn.ln()
    z = sum(D(j)**(-s) for j in range(1,n))
    dz = -sum(D(j).ln()*D(j)**(-s) for j in range(1,n))
    t = nn**(1-s)/(ss-1)
    z += t+nn**(-s)/2
    dz += t*(-ln-1/(ss-1))-ln*nn**(-s)/2
    for k in range(1,17):
        r = math.prod(range(s,s+2*k-1))
        t = dec(B[2*k]/math.factorial(2*k))*r*nn**(-s-2*k+1)
        z += t
        dz += t*(sum(1/D(s+j) for j in range(2*k-1))-ln)
    return z,dz

def arch(s):
    if s%2 == 0:
        psi = -GAMMA+sum(1/D(j) for j in range(1,s//2))
    else:
        psi = -GAMMA-2*D(2).ln()+2*sum(1/D(2*j-1) for j in range(1,(s-1)//2+1))
    return 1/D(s)+1/D(s-1)-PI.ln()/2+psi/2

def coeff(n):
    return [(-1)**(n+j)*math.comb(n,j)*math.comb(n+j,j) for j in range(n+1)]

def compress(values, degrees):
    cs = [coeff(n) for n in degrees]
    return [[sum(D(ci*cj)*(values[i]+values[j])/D(i+j+3)
                 for i,ci in enumerate(a) for j,cj in enumerate(b))
             for b in cs] for a in cs]

def pivots(a):
    n = len(a)
    l = [[D(0)]*n for _ in range(n)]
    ds = []
    for i in range(n):
        for j in range(i):
            l[i][j]=(a[i][j]-sum(l[i][k]*l[j][k]*ds[k] for k in range(j)))/ds[j]
        l[i][i]=1
        ds.append(a[i][i]-sum(l[i][k]**2*ds[k] for k in range(i)))
    return ds

def run():
    with localcontext() as ctx:
        ctx.prec=60
        degrees=list(range(2,11))
        for n in degrees:
            c=coeff(n)
            assert sum(F(x,j+1) for j,x in enumerate(c)) == 0
            assert sum(F(x,j+2) for j,x in enumerate(c)) == 0
        vs=[]
        for cutoff in [64,128]:
            vals=[]
            for s in range(2,13):
                z,dz=zeta_pair(s,cutoff)
                vals.append(arch(s)+dz/z)
            vs.append(vals)
        total=compress(vs[1],degrees)
        archmat=compress([arch(s) for s in range(2,13)],degrees)
        pm=[[archmat[i][j]-total[i][j] for j in range(len(degrees))] for i in range(len(degrees))]
        result={'status':'exact rational pole cancellation; Decimal diagnostics, not interval/Lean certificates',
                'degrees':degrees,'precision':60,
                'em_64_vs_128_max_difference':str(max(abs(a-b) for a,b in zip(*vs))),
                'total_LDL_pivots':[str(x) for x in pivots(total)],
                'arch_LDL_pivots':[str(x) for x in pivots(archmat)],
                'prime_LDL_pivots':[str(x) for x in pivots(pm)],
                'first_vector':coeff(2),'first_energy':str(total[0][0])}
        Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))

if __name__ == '__main__':
    run()
