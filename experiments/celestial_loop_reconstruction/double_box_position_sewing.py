"""Full massive Euclidean double box via exact position-space sewing.
All seven propagators retained. Not a celestial cut reconstruction.
"""
import numpy as np
from scipy.special import roots_legendre, roots_laguerre, kve, j1
from scipy.integrate import quad
from scipy.special import ndtri
from scipy.stats import qmc
import json

P1=np.array([.4,.2,0,0]); P2=np.array([-.1,.3,.2,0])
P4=np.array([.2,-.2,.1,.3]); P=P1+P2

def triangles(n,p,q):
    t,w=roots_legendre(n); t=(t+1)/2; w=w/2
    u,v=np.meshgrid(t,t,indexing="ij")
    wu,wv=np.meshgrid(w,w,indexing="ij")
    x=u.ravel(); y=((1-u)*v).ravel()
    wt=(wu*wv*(1-u)).ravel()
    shift=x[:,None]*p+y[:,None]*q
    mass=np.sqrt(1+x*np.dot(p,p)+y*np.dot(q,q)-np.sum(shift*shift,axis=1))
    return shift,mass,wt

def radial(ma,mb,q,n):
    t,w=roots_laguerre(n)
    scale=1+ma+mb
    r=t[None,:]/scale[:,None]
    z=q[:,None]*r
    # J1(qr)/q = r*(J1(z)/z); removable q=0 limit.
    ratio=np.full_like(z,.5)
    np.divide(j1(z),z,out=ratio,where=z!=0)
    value=r**4*ratio*kve(1,r)*kve(1,ma[:,None]*r)*kve(1,mb[:,None]*r)
    return (value@w)/(scale*ma*mb)

def sewn(n,nr):
    a,ma,wa=triangles(n,P1,P)
    b,mb,wb=triangles(n,-P4,P)
    total=0.
    for i in range(len(ma)):
        q=np.linalg.norm(b-a[i],axis=1)
        total+=wa[i]*np.dot(wb,radial(np.full_like(mb,ma[i]),mb,q,nr))
    return total

def direct_radial(ma,mb,q):
    def fun(r):
        if r==0:return 0.
        ratio=.5 if q==0 else j1(q*r)/(q*r)
        return r**4*ratio*kve(1,r)*kve(1,ma*r)*kve(1,mb*r)*np.exp(-(1+ma+mb)*r)/(ma*mb)
    return quad(fun,0,np.inf,epsabs=2e-12,epsrel=2e-12)

def parameter_check(power=18):
    means=[]
    for seed in range(8):
        z=qmc.Sobol(7,scramble=True,seed=8290+seed).random_base2(power)
        y=ndtri(z)**2; x=y/y.sum(axis=1)[:,None]
        h=x[:,6]; a=x[:,:3].sum(1)+h; b=x[:,3:6].sum(1)+h
        U=a*b-h*h
        v=x[:,1,None]*P1+x[:,2,None]*P
        w=-x[:,4,None]*P4+x[:,5,None]*P
        C=1+x[:,1]*np.dot(P1,P1)+x[:,2]*np.dot(P,P)+x[:,4]*np.dot(P4,P4)+x[:,5]*np.dot(P,P)
        F=U*C-b*np.sum(v*v,1)-a*np.sum(w*w,1)-2*h*np.sum(v*w,1)
        assert np.all(F>0)
        vals=(16*np.pi**3/15)*np.sqrt(np.prod(x,axis=1))*U/F**3
        means.append(float(vals.mean()))
    return dict(mean=float(np.mean(means)),empirical_se=float(np.std(means,ddof=1)/np.sqrt(8)),
                scrambles=8,power=power,values=means)

if __name__=="__main__":
    checks=[]
    for ma,mb,q in ((1.,1.,0.),(1.1,1.2,.4),(1.03,1.07,.8)):
        direct,err=direct_radial(ma,mb,q)
        approx=float(radial(np.array([ma]),np.array([mb]),np.array([q]),192)[0])
        assert abs(direct-approx)<2e-8
        checks.append(dict(ma=ma,mb=mb,q=q,adaptive=direct,error_estimate=err,laguerre192=approx))
    values=[]
    for n,nr in ((6,48),(10,96),(14,192)):
        val=sewn(n,nr)
        values.append(dict(triangle_order=n,radial_order=nr,value=val))
        print(json.dumps(values[-1]),flush=True)
    assert abs(values[-1]["value"]-.03671283863657149)<5e-6
    varied=[]
    originals=(P1.copy(),P2.copy(),P4.copy())
    for scale in (.5,1.,2.):
        P1,P2,P4=(scale*x for x in originals); P=P1+P2
        val=sewn(10,96); ref=parameter_check()
        assert abs(val-ref["mean"])<5*ref["empirical_se"]+2e-6
        row=dict(momentum_scale=scale,sewn=val,parameter=ref)
        varied.append(row); print(json.dumps(row),flush=True)
    with open("double_box_position_sewing_results.json","w") as f:
        json.dump(dict(radial_checks=checks,full_integral=values,varied_kinematics=varied,
            prior_parameter_estimate=.03671283863657149,
            prior_parameter_empirical_se=7.194839103379658e-7,
            scope="all-massive off-shell Euclidean scalar double box; numerical, not certified"),f,indent=2)
