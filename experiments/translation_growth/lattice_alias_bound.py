"""Exact-formula alias lower bounds for a fixed lattice of theta translates.

Two evaluations separated by 2*pi/step agree on every lattice translation.
Their xi-weighted difference annihilates the entire lattice trial space.
Evaluations use 40-digit mpmath; outputs are not interval certificates.
"""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps=40
TAU=mp.mpf('0.01')


def xi(s):
    return s*(s-1)/2 * mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)


def moment(w):
    return mp.pi*w/mp.sin(mp.pi*w) if w else mp.mpf(1)


def bound(step,delta,offset,t):
    omega=2*mp.pi/step
    z1=-delta+1j*(-omega/2+offset)
    z2=z1+1j*omega
    a,b=xi(mp.mpf('.5')+z2),-xi(mp.mpf('.5')+z1)
    size=max(abs(a),abs(b));a/=size;b/=size
    norm2=abs(a)**2*moment(2*z1.real)+abs(b)**2*moment(2*z2.real)
    norm2+=2*mp.re(a*mp.conj(b)*moment(z1+mp.conj(z2)))
    val=a*mp.exp(TAU*z1*z1-z1*t)+b*mp.exp(TAU*z2*z2-z2*t)
    # Cancellation for lattice basis elements is algebraic; test j=-7,0,11.
    error=max(abs(a*xi(mp.mpf('.5')+z1)*mp.exp(z1*j*step)+
                  b*xi(mp.mpf('.5')+z2)*mp.exp(z2*j*step))/size
              for j in [-7,0,11])
    return abs(val)/mp.sqrt(norm2),error


if __name__=='__main__':
    rows=[]
    for denominator in [12,16]:
        step=mp.mpf(1)/denominator
        best=(mp.mpf(0),None)
        for di in range(1,50):
            delta=mp.mpf(di)/100
            for oi in range(-4,5):
                offset=mp.mpf(oi)/4
                lower,error=bound(step,delta,offset,24)
                if lower>best[0]:best=(lower,(delta,offset,error))
        lower,(delta,offset,error)=best
        row=dict(step=str(step),t=24,lower_bound=str(lower),delta=str(delta),
                 offset=str(offset),annihilation_relative_error=str(error))
        rows.append(row);print(json.dumps(row),flush=True)
    Path(__file__).with_name('lattice_alias_results.json').write_text(json.dumps(rows,indent=2)+'\n')
