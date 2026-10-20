"""Finite zero-list diagnostic of prime-assembled packet minimizers; not RH proof."""
import struct,json,hashlib
from pathlib import Path
import mpmath as mp
from packet_cross import cross
from modulated_interval import events
mp.mp.dps=85
raw=Path(__file__).with_name('zeros_14.dat').read_bytes();pos=8
blocks=struct.unpack_from('<Q',raw)[0];zeros=[];counts=[]
for k in range(blocks):
 t0,t1,n0,n1=struct.unpack_from('<ddQQ',raw,pos);pos+=32;z=0
 for j in range(n1-n0):
  z+=int.from_bytes(raw[pos:pos+13],'little');pos+=13
  zeros.append(mp.mpf(t0)+mp.mpf(z)*mp.mpf(2)**-101)
 counts.append([t0,t1,n0,n1])
assert pos==len(raw)
assert all(a<b for a,b in zip(zeros,zeros[1:]))
checks=[mp.nstr(abs(zeros[j]-mp.im(mp.zetazero(j+1))),8) for j in range(3)]
rows=[]
for Ls in ['1.4','2','4']:
 L=mp.mpf(Ls);ev=events(L);ns=list(range(-8,9));ts=[j*mp.pi/L for j in ns]
 Q=mp.matrix([[cross(L,t,s,ev) for s in ts] for t in ts]);vals,V=mp.eigsy(Q);v=V[:,0]
 def ft(g):
  return mp.sqrt(2*L)*mp.fsum(v[k]*mp.sin((ts[k]+g)*L)/((ts[k]+g)*L) for k in range(17))
 terms=[abs(ft(g))**2+abs(ft(-g))**2 for g in zeros]
 total=mp.fsum(terms)
 endpoint=mp.fsum((-1)**j*v[k] for k,j in enumerate(ns))
 H=mp.mpf(5000)
 tail_est=endpoint**2*(mp.log(H/(2*mp.pi))+1)/(mp.pi*L*H)
 # tail bound per zero using the actual endpoint cancellation in the rational Fourier formula
 # no sum over omitted zeros is estimated or certified here
 rows.append({'L':Ls,'prime_min':mp.nstr(vals[0],24),'listed_zero_sum':mp.nstr(total,24),'fraction':mp.nstr(total/vals[0],15),'partial_fractions':{str(k):mp.nstr(mp.fsum(terms[:k])/vals[0],15) for k in [8,20,100,1000]},'heuristic_tail_ratio':mp.nstr(tail_est/(vals[0]-total),15),'parity_error':mp.nstr(mp.norm(v-mp.matrix(list(reversed(list(v))))),8),'first_eight_transform_abs':[mp.nstr(abs(ft(g)),10) for g in zeros[:8]]})
print(json.dumps({'scope':'Prime-side minimizer vs 4520 listed on-line zeros; non-certified diagnostic; omitted zeros unbounded here','sha256':hashlib.sha256(raw).hexdigest(),'blocks':counts,'first_three_height_errors':checks,'rows':rows},indent=2))
