"""Common-scale Mellin transform of the full massive off-shell double box.
Not the complete celestial amplitude; external ratios/angles remain fixed.
"""
import json
import numpy as np
from scipy.special import ndtri, gamma, roots_legendre
from scipy.stats import qmc
import double_box_position_sewing as pos

BASE=(pos.P1.copy(),pos.P2.copy(),pos.P4.copy())
S=3+.7j

def independent_position_mellin(order,triangle_order=6):
    z,w=roots_legendre(order); z=(z+1)/2; w=w/2
    t=np.tan(np.pi*z/2); jac=np.pi*(1+t*t)/2
    vals=[]
    for v in t:
        pos.P1,pos.P2,pos.P4=(v*p for p in BASE)
        pos.P=pos.P1+pos.P2
        vals.append(pos.sewn(triangle_order,96))
    return np.dot(w*jac*t**(S-1),vals)

def parameters(power=18):
    p1,p2,p4=BASE; P=p1+p2
    transformed=[]; masses=[1.,.1,.01,.001,0.]
    mass_values=[]
    for seed in range(8):
        z=qmc.Sobol(7,scramble=True,seed=9200+seed).random_base2(power)
        y=ndtri(z)**2; x=y/y.sum(1)[:,None]
        h=x[:,6]; a=x[:,:3].sum(1)+h; b=x[:,3:6].sum(1)+h
        U=a*b-h*h
        v=x[:,1,None]*p1+x[:,2,None]*P
        w=-x[:,4,None]*p4+x[:,5,None]*P
        C0=x[:,1]*np.dot(p1,p1)+x[:,2]*np.dot(P,P)+x[:,4]*np.dot(p4,p4)+x[:,5]*np.dot(P,P)
        F0=U*C0-b*np.sum(v*v,1)-a*np.sum(w*w,1)-2*h*np.sum(v*w,1)
        assert np.all(F0>0)
        measure=(8*np.pi**3/15)*np.sqrt(np.prod(x,axis=1))
        mellin=gamma(S/2)*gamma(3-S/2)/2*U**(S/2-2)*F0**(-S/2)
        transformed.append(complex(np.mean(measure*mellin)))
        row=[float(np.mean(measure*2*U/(m*U+F0)**3)) for m in masses]
        assert all(row[i]<=row[i+1] for i in range(len(row)-1))
        mass_values.append(row)
    def summary(v):
        v=np.asarray(v)
        return dict(mean=float(v.mean()),empirical_se=float(v.std(ddof=1)/np.sqrt(len(v))))
    return dict(power=power,scrambles=8,
        mellin_real=summary(np.real(transformed)),mellin_imag=summary(np.imag(transformed)),
        mass_scan=[dict(mass_squared=m,**summary(np.array(mass_values)[:,i])) for i,m in enumerate(masses)])

if __name__=="__main__":
    out=dict(s_real=S.real,s_imag=S.imag,parameter=parameters(),position=[])
    print(json.dumps(out["parameter"]),flush=True)
    out['initial_failed_refinement']=[dict(scale_nodes=24,triangle_order=6,real=.26430111151869806,imag=.4069140622108953),dict(scale_nodes=48,triangle_order=6,real=.2643020017890966,imag=.4069125231156313),dict(scale_nodes=96,triangle_order=6,real=.2643020309930757,imag=.40691252024938596)]
    for order,nt in ((32,10),(48,16)):
        value=independent_position_mellin(order,nt)
        row=dict(scale_nodes=order,triangle_order=nt,radial_order=96,real=value.real,imag=value.imag)
        out["position"].append(row)
        print(json.dumps(row),flush=True)
        with open('double_box_scale_mellin_results.json','w') as f: json.dump(out,f,indent=2)
    for key,part in (("mellin_real","real"),("mellin_imag","imag")):
        ref=out["parameter"][key]
        assert abs(out["position"][-1][part]-ref["mean"])<6*ref["empirical_se"]+2e-5
    with open("double_box_scale_mellin_results.json","w") as f: json.dump(out,f,indent=2)
