"""UNTRUSTED floating design. No continuum certificate from this script."""
import json, numpy as np
from scipy.special import digamma
from scipy.linalg import solve
import sys
D=json.load(open('even_input.json'));h=float(D['h']);C=D['C_num']/D['den'];a=np.array(D['a_num'])/D['den'];u=2*h+np.arange(len(a))*D['D'];dd=D['Dd'];kap=np.log(2)/np.sqrt(2)
N=int(sys.argv[1]) if len(sys.argv)>1 else 100
H=float(sys.argv[3]) if len(sys.argv)>3 else 6.;Xi=1500.;dx=float(sys.argv[2]) if len(sys.argv)>2 else .025
z=h+(H-h)*np.linspace(0,1,N+1)**1.5
# Unknowns lambda and N exterior nodal values; terminal value zero.
# Parseval surrogate replaces 1/m by 1/mX above Xi; mX below is a candidate analytic tail bound.
mX=float(digamma(.25+.5j*Xi).real-np.log(np.pi)-2*kap-3.3-4*C/(dd*Xi*Xi))
assert mX>0
M=np.zeros((N+1,N+1));v=np.zeros(N+1);s=0.
for start in np.arange(dx/2,Xi,2000*dx):
 x=np.arange(start,min(start+2000*dx,Xi),dx)
 mv=digamma(.25+.5j*x).real-np.log(np.pi)-2*kap*np.cos(x*np.log(2))+C*2*(2*h+dd/2)*np.sinc((2*h+dd/2)*x/np.pi)*np.sinc(dd*x/2/np.pi)-np.cos(x[:,None]*u)@a
 assert mv.min()>0
 base=np.sin(h*x)/x
 cw=(.5*np.sinh(h/2)*np.cos(h*x)+x*np.cosh(h/2)*np.sin(h*x))/(.25+x*x)
 co=np.cos(x[:,None]*z)
 Phi=np.empty((len(x),N)); dz=np.diff(z)
 Phi[:,0]=(co[:,0]-co[:,1])/(dz[0]*x*x)-np.sin(h*x)/x
 for j in range(1,N):Phi[:,j]=((co[:,j]-co[:,j-1])/dz[j-1]+(co[:,j]-co[:,j+1])/dz[j])/(x*x)
 R=np.column_stack([-cw,Phi]);weights=4/np.pi*(1/mv-1/mX)*dx
 M+=R.T@(weights[:,None]*R);v+=R.T@(weights*base);s+=base@(weights*base)
# Full Parseval mass for even residual extension, beta =2/mX int_0^inf r^2 plus correction
s+=2*h/mX;v[0]+=-4*np.sinh(h/2)/mX
M[0,0]+=(h+np.sinh(h))/mX
for j,dzj in enumerate(np.diff(z)):
 M[j+1,j+1]+=2*dzj/(3*mX)
 if j+1<N:
  M[j+2,j+2]+=2*dzj/(3*mX)
  M[j+1,j+2]+=dzj/(3*mX);M[j+2,j+1]+=dzj/(3*mX)
M[0,0]+=.5 # pole coordinate cost
sol=solve(M,-v,assume_a='pos');val=s+v@sol
print(json.dumps(dict(N=N,dx=dx,lambda_value=sol[0],cost=val,C=C,margin=1-C*val,mX=mX,min_eig=np.linalg.eigvalsh(M)[0]),indent=2),flush=True)
json.dump(dict(status='UNTRUSTED FLOAT DESIGN ONLY',h=h,H=H,Xi=Xi,dx=dx,knots=z.tolist(),lambda_value=sol[0],values=sol[1:].tolist(),cost=val,margin=1-C*val),open(f'even_design_N{N}.json','w'))
