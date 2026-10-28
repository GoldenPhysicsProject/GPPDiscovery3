"""Classical polarized Rosen vacuum wave: nonlinear reconstruction, not loop amplitudes."""
import json, math
import numpy as np
from scipy.integrate import solve_ivp

# Coordinates (u,v,x,y), signature -+++; initial data a(-4)=1,a'(-4)=0.
def profile(u):
    f=np.exp(-u*u)
    return f,-2*u*f,(4*u*u-2)*f

def ricci(a,ap,app,f,fp,fpp,e):
    g=np.zeros((4,4)); g[0,1]=g[1,0]=-1
    gp=np.zeros_like(g); gpp=np.zeros_like(g)
    for i,sgn in [(2,1),(3,-1)]:
        ex=np.exp(sgn*e*f)
        g[i,i]=a*a*ex
        gp[i,i]=ex*(2*a*ap+sgn*e*fp*a*a)
        gpp[i,i]=ex*(2*ap*ap+2*a*app+4*sgn*e*fp*a*ap+sgn*e*fpp*a*a+e*e*fp*fp*a*a)
    inv=np.linalg.inv(g); invp=-inv@gp@inv
    dg=np.zeros((4,4,4)); ddg=np.zeros_like(dg)
    dg[0]=gp; ddg[0]=gpp
    G=np.zeros((4,4,4)); Gp=np.zeros_like(G)
    for r in range(4):
      for m in range(4):
       for n in range(4):
        for s in range(4):
         b=dg[m,s,n]+dg[n,s,m]-dg[s,m,n]
         bp=ddg[m,s,n]+ddg[n,s,m]-ddg[s,m,n]
         G[r,m,n]+=.5*inv[r,s]*b
         Gp[r,m,n]+=.5*(invp[r,s]*b+inv[r,s]*bp)
    R=np.zeros((4,4))
    for m in range(4):
     for n in range(4):
      R[m,n]=Gp[0,m,n]-(sum(Gp[r,m,r] for r in range(4)) if n==0 else 0)
      for r in range(4):
       for s in range(4):
        R[m,n]+=G[r,m,n]*G[s,r,s]-G[s,m,r]*G[r,n,s]
    return R

def run(e):
    N=20
    def hierarchy(u,y):
        dy=np.zeros_like(y); dy[::2]=y[1::2]
        q=profile(u)[1]**2/4
        dy[3::2]=-q*y[:-2:2]
        return dy
    y=np.zeros(2*(N+1)); y[0]=1
    hs=solve_ivp(hierarchy,(-4,2),y,rtol=2e-12,atol=1e-14,dense_output=True,method='DOP853')
    def ode(u,y):return [y[1],-e*e*profile(u)[1]**2*y[0]/4]
    sol=solve_ivp(ode,(-4,2),[1,0],rtol=2e-12,atol=1e-14,dense_output=True,method='DOP853')
    us=np.linspace(-4,2,161); exact=sol.sol(us)[0]
    assert exact.min()>.1, 'Rosen caustic: coordinate chart invalid'
    cs=hs.sol(us)[::2]; rows=[]
    for n in [1,2,5,10,20]:
        estimate=sum(e**(2*k)*cs[k] for k in range(n+1))
        # sup f'^2 = 2/euler; q=f'^2/4 <=1/(2euler)
        B=e*e/(2*math.e); T=6
        # Bound tail of cosh(sqrt(B)*T) by positive terms, with geometric remainder.
        first=(B*T*T)**(n+1)/math.factorial(2*n+2)
        ratio=B*T*T/((2*n+3)*(2*n+4))
        assert ratio<1
        bound=first/(1-ratio)
        rows.append(dict(order_in_amplitude=2*n,max_error=float(np.max(abs(estimate-exact))),analytic_series_tail_bound=bound))
    residual=0.; bad=0.
    for u in us[::8]:
        a,ap=sol.sol(u); f,fp,fpp=profile(u); app=-e*e*fp*fp*a/4
        residual=max(residual,float(np.max(abs(ricci(a,ap,app,f,fp,fpp,e)))))
        bad=max(bad,abs(ricci(1,0,0,f,fp,fpp,e)[0,0]))
    assert residual<1e-10 and rows[-1]['max_error']<1e-9 and bad>1e-3
    return dict(amplitude=e,minimum_a=float(exact.min()),ricci_residual=residual,uncorrected_ricci=bad,series=rows)
if __name__=='__main__':
    results=[run(e) for e in [.3,.6,1.]]
    print(json.dumps(results,indent=2))
    open('rosen_nonlinear_reconstruction_results.json','w').write(json.dumps(results,indent=2)+'\n')
