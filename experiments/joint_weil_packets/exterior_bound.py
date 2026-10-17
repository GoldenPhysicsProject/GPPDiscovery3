"""Explicit Fourier exterior bound budget. Numerical evaluation of analytic formulas;
not interval certification, not RH, and no enormous matrix is allocated.
"""
import mpmath as mp,json
from modulated_interval import events
mp.mp.dps=50
a0=mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi)
rows=[]
for Ls in ['0.5','1.4','2','3','4']:
    L=mp.mpf(Ls)
    # Recover sum of the two translation weights from triangle-weighted events.
    s=mp.fsum(w/(1-u/(2*L)) for u,w in events(L))
    b=s+2; T=2*mp.pi*mp.exp(s+3)
    a=L*T/mp.pi
    cutoff=max(a+4*L*T*(b-a0)/mp.pi**3,16*L*mp.sinh(L/2)**2/mp.pi**2,mp.mpf(1))
    N=mp.floor(cutoff)+1
    leakage=4*L*T/(mp.pi**3*(N-a))
    negative_pole=8*L*mp.sinh(L/2)**2/(mp.pi**2*N)
    lower=b-(b-a0)*leakage-s-negative_pole
    assert leakage<=1/(b-a0) and negative_pole<=mp.mpf('.5')
    # Rounding evaluation only; proof of these inequalities is analytic.
    assert lower>=mp.mpf('.49999999999999999999')
    rows.append(dict(L=Ls,s=mp.nstr(s,18),log10_T=mp.nstr(mp.log10(T),18),log10_N=mp.nstr(mp.log10(N),18),lower=mp.nstr(lower,18)))
# Independent finite diagnostic of the actual Fourier low-band leakage trace.
controls=[]
for L,T,N in [(mp.mpf(1),mp.mpf(2),4),(mp.mpf('1.4'),mp.mpf(3),7)]:
    finite=mp.quad(lambda t:mp.fsum(2*L*mp.sin(L*t)**2/((L*t-j*mp.pi)**2) for j in list(range(N+1,81))+list(range(-80,-N)))/(2*mp.pi),[-T,0,T])
    bound=4*L*T/(mp.pi**3*(N-L*T/mp.pi))
    assert finite<bound
    controls.append(dict(L=str(L),T=str(T),N=N,partial_trace=mp.nstr(finite,18),analytic_upper=mp.nstr(bound,18)))
print(json.dumps(dict(scope='Explicit positive exterior bound, not finite-block positivity. Numerical budget diagnostics only.',rows=rows,controls=controls),indent=2))
