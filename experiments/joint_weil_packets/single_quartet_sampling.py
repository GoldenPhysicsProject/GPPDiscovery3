"""Synthetic quartet injection, not an alteration of zeta or a prime-side proof."""
import json
import mpmath as m
m.mp.dps=60

def B(z):
    return 8*m.pi**4*(m.cosh(z)-1)/(z**2*(z**2+4*m.pi**2)**2)
def quartet(u,z,b):
    return 4*m.re(b*m.cosh(z*u))
rows=[]
M=m.mpf('.001')
for delta,gamma in [('.1','14.134725'),('.1','100'),('.01','100'),('.1','1000')]:
    z=m.mpc(delta,gamma); b=B(z); c0=M+4*m.re(b)
    # Select positive phases of the exponentially growing term.
    bound=m.log((abs(c0)+M+2*abs(b)+1)/(2*abs(b)))/m.re(z)
    n=int(m.ceil((m.im(z)*bound+m.arg(b))/(2*m.pi)))
    u=(2*m.pi*n-m.arg(b))/m.im(z)
    q=quartet(u,z,b)
    lower=2*abs(b)*(m.exp(m.re(z)*u)-m.exp(-m.re(z)*u))
    assert q>=lower-m.mpf('1e-45')*max(1,abs(q))
    cu=M+q
    mineig=c0-abs(cu)
    assert mineig<0
    rows.append(dict(delta=delta,gamma=gamma,B_abs=str(abs(b)),lag=str(u),negative_eigenvalue=str(mineig)))
# Deliberate exact blind-grid model with b=i (not claiming b=B(z)).
h=m.mpf('.7');z=m.mpc('.1',m.pi/h);b=1j
blind=[quartet(n*h,z,b) for n in range(21)]
assert max(map(abs,blind))<m.mpf('1e-55')
other=[quartet(n*h*m.sqrt(2),z,b) for n in range(21)]
assert max(map(abs,other))>1
print(json.dumps({"scope": __doc__, "actual_packet_multiplier_examples": rows,
 "alias_control": {"first_grid_max": str(max(map(abs,blind))),
                   "second_grid_max": str(max(map(abs,other)))}}, indent=2))
