"""Diagnostics of joint-transform cancellation and sharp Abel pole bound.
Not interval arithmetic, not a zeta-zero certificate, not an RH proof.
"""
import json
import mpmath as mp
mp.mp.dps=70
rows=[]
for z in [mp.mpc('0.8','0.7'),mp.mpc('1.2','2.3'),mp.mpc('2','0')]:
 s=z+mp.mpf('.5'); Z=mp.zeta(s); Zd=mp.diff(mp.zeta,s)
 D=-Zd/(z*s)-2/(z-mp.mpf('.5'))**2+6/(z-mp.mpf('.5'))
 joint=(2/(z-mp.mpf('.5'))**2-6/(z-mp.mpf('.5'))+D)/Z-2/(z-mp.mpf('.5'))+2/z
 direct=-(Zd/Z+1/(z-mp.mpf('.5'))+1)/(z*s)
 err=abs(joint-direct)
 assert err<mp.mpf('1e-60')
 rows.append({'z':str(z),'transform_error':str(err)})
# f(t)=2 Re(r exp((a0+i gamma)t)); exact Abel norm and transform.
r=mp.mpc('0.7','-0.2'); gamma=mp.mpf('3'); a0=mp.mpf('.13')
energy=[]
for h in map(mp.mpf,['.1','.01','.001','.0001','.00001']):
 A=abs(r)**2/h+mp.re(r*r/(h-1j*gamma))
 F=r/(2*h)+mp.conj(r)/(2*h+2j*gamma)
 lower=4*h*abs(F)**2/(1+h/mp.sqrt(h*h+gamma*gamma))
 assert 0<=lower<=A+mp.mpf('1e-60')
 energy.append({'h':str(h),'h_A':str(h*A),'h_lower':str(h*lower),'limit_r_squared':str(abs(r)**2)})
print(json.dumps({'scope':__doc__,'joint_transform':rows,'synthetic_sharpness':energy},indent=2))
