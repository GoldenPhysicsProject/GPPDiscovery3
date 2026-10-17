"""Checks the analytic 1/|j-i| joint-form envelope. No interval certification."""
import mpmath as mp,json
from packet_cross import cross
from modulated_interval import events
mp.mp.dps=65
rows=[]
for Ls in ['0.5','1.4','2']:
    L=mp.mpf(Ls); ev=events(L); N=3
    s=mp.fsum(w/(1-u/(2*L)) for u,w in ev)
    H=2+mp.pi/2+4*mp.exp(-L)/(1-mp.exp(-4*L))
    worst=mp.mpf(0)
    for i in range(-N,N+1):
        ci=(H+s)/mp.pi+4*mp.sinh(L/2)**2/(mp.pi*mp.sqrt(mp.mpf(1)/4+(i*mp.pi/L)**2))
        for j in [-12,-7,-4,4,7,12]:
            value=abs(cross(L,j*mp.pi/L,i*mp.pi/L,ev))
            bound=ci/(abs(j)-N)
            assert value<bound
            worst=max(worst,value/bound)
    rows.append(dict(L=Ls,tests=42,max_ratio=mp.nstr(worst,20)))
print(json.dumps(dict(scope='Numerical checks of analytic coupling envelope; proof is separate.',rows=rows),indent=2))
