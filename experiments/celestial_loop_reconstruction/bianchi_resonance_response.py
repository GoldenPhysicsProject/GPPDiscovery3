import mpmath as m,json
m.mp.dps=45
beta=lambda s:m.dirichlet(s,[0,1,0,-1])
L=lambda s:m.power(m.pi,-s)*m.gamma(s)*m.zeta(s)*beta(s)
Phi=lambda s:L(s-1)/L(s)
rho=m.zetazero(1); gamma=m.im(rho)
rows=[]; maxunit=m.mpf(0); maxder=m.mpf(0)
for offset in [-1,-.5,0,.5,1]:
 t=gamma+offset; s=1+1j*t
 ph=Phi(s)
 response=2*m.re(m.diff(L,s)/L(s))
 direct=m.re(1j*m.diff(lambda x:Phi(1+1j*x),t)/ph)
 maxunit=max(maxunit,abs(abs(ph)-1)); maxder=max(maxder,abs(response-direct))
 rows.append({'lambda':float(t),'oriented_full_phase_response':float(response),'isolated_rho_factor':float(1/(m.mpf('.25')+offset**2))})
residue=L(rho-1)/m.diff(L,rho)
reschecks=[]
for eps in [m.mpf('1e-4'),m.mpf('1e-6'),m.mpf('1e-8')]:
 err=abs(eps*Phi(rho+eps)/residue-1)
 reschecks.append({'epsilon':float(eps),'relative_residue_error':float(err)})
assert maxunit<m.mpf('1e-35') and maxder<m.mpf('1e-35')
assert reschecks[-1]['relative_residue_error']<1e-6
# Pole/zero factor with arbitrary off-line beta: positivity survives.
controls=[]
for b in [m.mpf('.3'),m.mpf('.5'),m.mpf('.7')]:
 a=1-b
 pole=b+1j*gamma
 B=lambda s:(s-(2-m.conj(pole)))/(s-pole)
 t=gamma+m.mpf('.2')
 val=m.re(1j*m.diff(lambda x:B(1+1j*x),t)/B(1+1j*t))
 formula=2*a/(a*a+(t-gamma)**2)
 assert abs(val-formula)<m.mpf('1e-35')
 controls.append({'beta':float(b),'response':float(val),'positive':bool(val>0)})
result={'precision_digits':45,'first_zeta_ordinate':float(gamma),'max_unitarity_error':float(maxunit),'max_derivative_identity_error':float(maxder),'samples':rows,'simple_pole_residue':{'real':float(m.re(residue)),'imag':float(m.im(residue))},'residue_checks':reschecks,'off_line_factor_controls':controls,'scope':'Floating checks of established scattering formula and derived local response; no rigorous interval certification or Lean verification.'}
print(json.dumps(result,indent=2))
