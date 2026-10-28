import numpy as np, mpmath as mp,json
from scipy.integrate import simpson
mp.mp.dps=25
gamma=float(mp.im(mp.zetazero(1)))
k=np.linspace(gamma-2,gamma+2,513); dk=k[1]-k[0]
a=np.exp(-.5*((k-gamma)/.25)**2); a/=np.sqrt(simpson(a*a,x=k))
w=np.sqrt(1+k*k)
L=lambda s:mp.pi**(-s)*mp.gamma(s)*mp.zeta(s)*mp.dirichlet(s,[0,1,0,-1])
# Functional equation: Phi=conj(A)/A; outgoing/incoming S=1/Phi.
S=np.array([complex(L(1+1j*float(v))/mp.conj(L(1+1j*float(v)))) for v in k])
phase=np.unwrap(np.angle(S)); phaseprime=np.gradient(phase,dk,edge_order=2)
phaseprime[2:-2]=(phase[:-4]-8*phase[1:-3]+8*phase[3:-1]-phase[4:])/(12*dk)
avgshift=float(simpson(a*a*phaseprime,x=k))
vmean=float(simpson(a*a*k/w,x=k))
r=np.linspace(-80,80,6401)
weights=np.ones(len(k));weights[1:-1:2]=4;weights[2:-1:2]=2;weights*=dk/3
kernel=np.exp(1j*r[:,None]*k[None,:])/np.sqrt(2*np.pi)
rows=[]
for t in (0.,20.,40.):
 coeff=a*S*np.exp(-1j*w*t)
 out=kernel@(weights*coeff)
 ref=kernel@(weights*a*np.exp(-1j*w*t))
 norm=float(simpson(abs(out)**2,x=r)); refnorm=float(simpson(abs(ref)**2,x=r))
 mean=float(simpson(r*abs(out)**2,x=r)/norm)
 refmean=float(simpson(r*abs(ref)**2,x=r)/refnorm)
 rows.append({'time':t,'out_norm':norm,'reference_norm':refnorm,'out_centroid':mean,'reference_centroid':refmean,'displacement':mean-refmean,'predicted_displacement':-avgshift})
 assert abs(norm-1)<1e-9 and abs(refnorm-1)<1e-9
 assert abs(mean-refmean+avgshift)<.002
# Independent finite-difference PDE check of a cusp-channel packet.
def u(t,x):return np.sum(weights*a*S*np.exp(1j*k*x-1j*w*t))/np.sqrt(2*np.pi)
h=1e-4; residuals=[]
for t,x in [(0,-5),(2,-3),(10,5)]:
 mid=u(t,x)
 utt=(u(t+h,x)-2*mid+u(t-h,x))/h**2
 urr=(u(t,x+h)-2*mid+u(t,x-h))/h**2
 residuals.append(float(abs(utt-urr+mid)))
assert max(residuals)<1e-5
result={'model':'ultrastatic R_time x PSL2(Z[i])\\H3, massless scalar; numerics are cusp-channel envelopes, not full interior Eisenstein evaluation','center':gamma,'sigma_amplitude':.25,'spectral_points':len(k),'max_unitarity_error':float(np.max(abs(abs(S)-1))),'average_group_velocity':vmean,'phase_shift_average':avgshift,'narrow_packet_delay_estimate':avgshift/vmean,'packets':rows,'finite_difference_wave_equation_max_residual':max(residuals),'milne_mode_ode_coefficient_max':float(max(abs((-1-1j*q)**2+2*(-1-1j*q)+1+q*q) for q in k)), 'integrated_flux_mismatch':float(abs(simpson(k*a*a*(abs(S)**2-1),x=k)))}
print(json.dumps(result,indent=2))
