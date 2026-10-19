"""Independent Fourier-packet audit of Claude's pole-free inertia observation.
Numerical, not an inertia certificate or RH proof.
"""
import json
import mpmath as mp
from packet_cross import cross
from modulated_interval import events
mp.mp.dps=85
rows=[]
for Ls in ['0.5','1.4','2','3','4']:
    L=mp.mpf(Ls);n=8;ev=events(L);ts=[j*mp.pi/L for j in range(-n,n+1)]
    Q=mp.matrix([[cross(L,t,s,ev) for s in ts] for t in ts])
    ell=lambda t:2*mp.sinh((mp.mpf(1)/2+1j*t)*L)/((mp.mpf(1)/2+1j*t)*mp.sqrt(2*L))
    P=mp.matrix([[2*mp.re(ell(t)*ell(s)) for s in ts] for t in ts])
    E=mp.matrix(2*n+1,n+1);E[n,0]=1
    O=mp.matrix(2*n+1,n)
    for j in range(1,n+1):
        E[n+j,j]=E[n-j,j]=1/mp.sqrt(2)
        O[n+j,j-1]=1/mp.sqrt(2);O[n-j,j-1]=-1/mp.sqrt(2)
    out={'L':Ls}
    for sec,V in [('even',E),('odd',O)]:
        R=V.T*(Q-P)*V;W=V.T*Q*V
        eigR=mp.eigsy(R,eigvals_only=True);eigW=mp.eigsy(W,eigvals_only=True)
        if sec=='even':
            v=mp.matrix([mp.sqrt(2)*mp.re(ell(0))]+[2*mp.re(ell(j*mp.pi/L)) for j in range(1,n+1)])
            gap=1+(v.T*mp.lu_solve(R,v))[0]
        else:
            v=mp.matrix([2*mp.im(ell(j*mp.pi/L)) for j in range(1,n+1)])
            gap=1-(v.T*mp.lu_solve(R,v))[0]
        sign=1 if sec=='even' else -1
        assert mp.norm(W-R-sign*v*v.T)<mp.mpf('1e-70')
        out[sec]=dict(dimension=R.rows,R_negative=sum(x<0 for x in eigR),Q_negative=sum(x<0 for x in eigW),min_R=mp.nstr(eigR[0],14),min_Q=mp.nstr(eigW[0],14),schur_margin=mp.nstr(gap,14),margin_required='<=0' if sec=='even' else '>=0')
    rows.append(out)
print(json.dumps({'scope':'Finite Fourier-packet subspaces, 85 decimal digits, analytic series tails not rounded by intervals.','rows':rows},indent=2))
