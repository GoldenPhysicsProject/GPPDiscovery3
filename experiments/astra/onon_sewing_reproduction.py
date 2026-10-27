"""Reproduction and repair of ONON shadow sewing.
Numerical tests, not Lean proofs or a complete multi-loop amplitude.
Run: python onon_sewing_reproduction.py
Dependencies: mpmath, scipy.
"""
import json, math
import mpmath as mp
from scipy.integrate import dblquad
mp.mp.dps = 40

def close(a,b,tol):
    assert abs(a-b) <= tol*max(1,abs(a),abs(b)), (a,b)

def sew_exponentials(a,b):
    a,b=mp.mpf(a),mp.mpf(b)
    def integrand(t):
        return (mp.gamma(1+1j*t)*mp.gamma(1-1j*t)
                *a**(-1-1j*t)*b**(-1+1j*t))
    spectral=mp.quad(integrand,[-mp.inf,0,mp.inf])/(2*mp.pi)
    energy=mp.quad(lambda w:w*mp.exp(-(a+b)*w),[0,mp.inf])
    exact=1/(a+b)**2
    close(spectral,exact,mp.mpf("1e-30"))
    close(energy,exact,mp.mpf("1e-30"))
    return {"a":str(a),"b":str(b),"spectral":str(spectral),
            "energy":str(energy),"absolute_error":str(abs(spectral-energy))}

def coupled_sewing():
    # f(w)=exp(-sum(log w))*exp(-x^T A x/2), and analogously B.
    # Mellin contours c=(1,1); non-diagonal A,B test genuine coupling.
    A=mp.matrix([[2,mp.mpf(".5")],[mp.mpf(".5"),1]])
    B=mp.matrix([[1,mp.mpf("-.25")],[mp.mpf("-.25"),3]])
    C=A+B; H=A**-1+B**-1
    assert mp.det(A)>0 and mp.det(B)>0
    # Independent 2D numerical quadratures, finite +/-12 limits.
    def quad_gauss(M):
        a,b,d=float(M[0,0]),float(M[0,1]),float(M[1,1])
        val,err=dblquad(lambda y,x:math.exp(-.5*(a*x*x+2*b*x*y+d*y*y)),
                        -12,12,-12,12,epsabs=1e-10,epsrel=1e-10)
        return val,err
    direct,de=quad_gauss(C)
    spectral,se=quad_gauss(H)
    pref=1/math.sqrt(float(mp.det(A)*mp.det(B)))
    spectral*=pref; se*=pref
    exact=2*mp.pi/mp.sqrt(mp.det(C))
    close(direct,float(exact),1e-9)
    close(spectral,float(exact),1e-9)
    return {"direct_log_energy":direct,"spectral":spectral,
            "exact":str(exact),"quadrature_error_estimates":[de,se],
            "finite_limits":[-12,12],"interval_certified":False}

def old_two_loop():
    def w(x):
        return 1.0 if x==0 else math.pi*x/math.sinh(math.pi*x)
    # Gamma(3/2)^4/Gamma(4)=pi^2/96; paired Gamma = w(x).
    for x,y in ((.2,.7),(-1.3,.4),(2,3)):
        gamma=mp.gamma(mp.mpf("1.5"))**4/mp.gamma(4)
        gamma*=mp.gamma(1+1j*x)*mp.gamma(1-1j*x)
        gamma*=mp.gamma(1+1j*y)*mp.gamma(1-1j*y)
        close(float(mp.re(gamma))*w(x)*w(y),
              math.pi**2/96*w(x)**2*w(y)**2,1e-13)
    pref=math.pi**2/96/(2*math.pi)**4
    results=[]
    for R in (4,8,12):
        val,err=dblquad(lambda y,x:w(x)**2*w(y)**2/(1+x*x+y*y),
                        -R,R,-R,R,epsabs=1e-10,epsrel=1e-10)
        results.append({"cutoff":R,"value":pref*val,"error_estimate":pref*err})
    # Integral w(x)^2 over R is pi/3. Drop denominator >=1 for exact bound.
    bound=1/13824
    assert results[-1]["value"] < bound < .0172267891
    return {"original_integrand_off_axes":results,
            "original_axis_bug":"returns 1 for both factors when either variable is near zero",
            "printed_claim":.0172267891,"printed_reference_decimal":.0172268049,
            "actual_reference_expression_pi6_over_5670":math.pi**6/5670,
            "corrected_continuous_integrand_exact_upper_bound":bound}

def G(c,d):
    w=d+4*c
    return 8*mp.atanh(mp.sqrt(d/w))/(mp.sqrt(d)*w**mp.mpf("1.5"))+2/(c*w)

def box_direct(s,u,m):
    S,U=-s,-u
    def inner(x,y):
        r=1-x-y
        if r<=0:return mp.mpf(0)
        return r*G(m*(x+y)+U*x*y,S*r*r)
    return mp.quad(lambda x:mp.quad(lambda v:(1-x)*inner(x,(1-x)*v),
                      [0,.5,1]),[0,.01,.5,.99,1])

def box_cut(v,u,m):
    c=m*(v+m);d=v*(-u);w=d+4*c
    return mp.atanh(mp.sqrt(d/w))/(2*mp.pi*mp.sqrt(d*w))

def box_disp(s,u,m):
    return 8*mp.pi*mp.quad(lambda v:box_cut(v,u,m)/(v-s),
                               [0,m,1,10,100,1000,mp.inf])

if __name__=="__main__":
    report={"sewing":[sew_exponentials(a,b) for a,b in
                     (("1","1"),("1.3","2.7"),(".2","4"))],
            "coupled_two_variable_sewing":coupled_sewing(),
            "old_two_loop_audit":old_two_loop(),"one_loop_box":[]}
    for point in (("-3","-2",".5"),("-2","-5",".7"),("-1.25","-.8",".2")):
        s,u,m=map(mp.mpf,point)
        d=box_direct(s,u,m);r=box_disp(s,u,m)
        close(d,r,mp.mpf("1e-12"))
        report["one_loop_box"].append({"s":str(s),"u":str(u),"mu2":str(m),
            "direct":str(d),"dispersion":str(r),"absolute_difference":str(abs(d-r))})
    print(json.dumps(report,indent=2))
