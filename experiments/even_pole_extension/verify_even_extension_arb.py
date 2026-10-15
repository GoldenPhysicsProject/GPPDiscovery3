"""Astra adaptation of Grok S25 Arb verifier for the EVEN pole-completion test.
Q=A+2|w><w|-C|1><1|; r=1-lambda*cosh(x/2) extended evenly.
S is its half-line COSINE transform, beta=4/pi integral S^2/m.
Certificate criterion C*(beta+lambda^2/2)<1. RH OPEN.
Original S25 cell, Bernstein and Parseval routines retained; independent review requested.
"""
import json, sys, time, os
from flint import arb, acb, acb_poly, fmpq, ctx
ctx.prec = int(os.environ.get("PREC", "128"))
REL = os.environ.get("REL", "0.005")
XBETA_EXACT = int(os.environ.get("XBE", "3000"))
LOG2R0 = 7  # initial half-width r0 = 2^-7 (dyadic, so all cell points are exact)

def make_S(h, ys, zs, lam):
    n=len(ys)-2
    sl=[(ys[i+1]-ys[i])/(zs[i+1]-zs[i]) for i in range(n+1)]
    coef=[(zs[0],-sl[0])]+[(zs[j],-(sl[j]-sl[j-1])) for j in range(1,n+1)]+[(zs[n+1],sl[n])]
    def S(x):
        if x == 0:
            return h-2*lam*(h/2).sinh()+sum(((ys[i]+ys[i+1])*(zs[i+1]-zs[i])/2 for i in range(n+1)),arb(0))
        cw=arb(0)
        for al in (arb(1)/2,-arb(1)/2):
            zz=acb(al,x);cw+=(((zz*h).exp()-1)/zz).real/2
        base=(h*x).sin()/x-lam*cw
        acc=sum((cf*(x*z).cos() for z,cf in coef),arb(0))
        return base-ys[0]*(x*h).sin()/x+acc/(x*x)
    return S,sl

def verify(fn):
    d = json.load(open(fn)); t0 = time.time(); ext = d["extension"]
    hq = fmpq(*d["h_q"]); h = arb(hq); Hq = fmpq(*ext["H_q"]); H = arb(Hq); N = ext["N"]
    A_ = [arb(fmpq(n, d["den"])) for n in d["a_num"]]; C = arb(fmpq(d["C_num"], d["den"]))
    c = fmpq(*d["c"]); Dsp = arb(fmpq(*d["D_q"])); Dd = arb(fmpq(*d["Dd_q"]))
    ys = [arb(fmpq(n, ext["alpha_den"])) for n in ext["alpha_num"]] + [arb(0)]
    if "z_num" in ext:
        zs = [arb(fmpq(n, ext["z_den"])) for n in ext["z_num"]]
        assert fmpq(ext["z_num"][0], ext["z_den"]) == hq and fmpq(ext["z_num"][-1], ext["z_den"]) == Hq
        assert all(ext["z_num"][i] < ext["z_num"][i+1] for i in range(N+1))
    else:  # S23/S24 uniform knots
        zs = [arb(hq + (Hq-hq)*fmpq(k, N+1)) for k in range(N+2)]
    assert len(zs) == len(ys) == N+2
    zs[0] = h; zs[-1] = H
    lam=arb(fmpq(*d["lambda_q"]))
    Sfun, sl = make_S(h, ys, zs, lam)
    # M3 >= int_0^H t^3 |v_w|,  ||v_w||^2 on (0, inf) exactly
    M3 = h**4/4*(abs(1-lam)+abs(lam)*((h/2).cosh()-1)) + sum(((zs[i+1]**3)*(abs(ys[i])+abs(ys[i+1]))*(zs[i+1]-zs[i])/2 for i in range(N+1)), arb(0))
    vnorm2 = (h-4*lam*(h/2).sinh()+lam*lam*(h+(h).sinh())/2) + sum(((zs[i+1]-zs[i])*(ys[i]**2+ys[i]*ys[i+1]+ys[i+1]**2)/3 for i in range(N+1)), arb(0))
    print(f"[info] N={N} M3={float(M3.mid()):.3f} ||v_w||^2={float(vnorm2.mid()):.6f} jump={float(abs(1-lam*(h/2).cosh()-ys[0]).mid()):.4f}", flush=True)
    P = acb_poly([acb(x) for x in A_]); dP = P.derivative()
    tlg = arb(2).log(); kap = tlg/arb(2).sqrt(); a0 = arb(1)/4; lpi = arb.pi().log()
    L0 = 2*h; u = [L0 + Dsp*k for k in range(d["K"])]; Aw = 2*h + Dd/2; bdc = Dd/2
    NN = 2**16; mx = arb(0)
    for j in range(NN//2+1):
        z = (acb(0, 2*arb.pi()*j/NN)).exp(); v = abs(P(z)).upper()
        if v > mx: mx = v
    ber = 1 - arb.pi()*(d["K"]-1)/NN; assert ber > 0; Mbern = (mx/ber).upper()
    Xi = arb(d["Xi_LP"]).floor() + 1
    def mlow_tail(X):
        return acb(a0, X/2).digamma().real - lpi - 2*kap - Mbern - C*4/(Dd*X*X)
    tail = mlow_tail(Xi) - arb(c); assert tail > 0
    z3 = acb(3).zeta(acb(a0)).real
    F2 = (z3/2 + 2*kap*tlg*tlg + C*2*(2*h+Dd)**3/3 + sum(abs(x)*uk*uk for x, uk in zip(A_, u))).upper()
    def fval(x):
        w = acb(a0, x/2); z = acb(0, Dsp*x).exp(); e = acb(0, L0*x).exp(); Pz = P(z)
        s1 = (Aw*x).sin()/x; s2 = (bdc*x).sin()/(bdc*x)
        s1p = (Aw*x*(Aw*x).cos() - (Aw*x).sin())/(x*x); s2p = (bdc*x*(bdc*x).cos() - (bdc*x).sin())/(bdc*x*x)
        fv = w.digamma().real - lpi - 2*kap*(x*tlg).cos() + C*2*s1*s2 - (e*Pz).real - arb(c)
        Sp = (acb(0, L0)*e*Pz + e*dP(z)*acb(0, Dsp)*z).real
        fp = -acb(2).zeta(w).imag/2 + 2*kap*tlg*(x*tlg).sin() + C*2*(s1p*s2 + s1*s2p) - Sp
        return fv, fp
    XBi = max(int(float(d["Xi_LP"])) + 2, XBETA_EXACT); XB = arb(XBi)
    nx = XBi*2**(LOG2R0-1); r0 = arb(1)/2**LOG2R0
    stack = [(r0*(2*i+1), r0) for i in range(nx)][::-1]
    cache = {}
    def Sc(x):
        k = int(round(float(x.mid())*2**34))
        if k not in cache: cache[k] = Sfun(x)
        return cache[k]
    ncell = nsplit = nref = 0; worst = None
    bint = arb(0); ssq_low = arb(0); rel = arb(REL); s3 = arb(3).sqrt()
    while stack:
        mid, r = stack.pop()
        fv, fp = fval(mid); slack = abs(fp)*r + F2*r*r/2; lb = fv - slack; ncell += 1
        if lb > 0 and slack > (fv+arb(c))*rel and r > arb("1e-5"):
            nref += 1; stack += [(mid+r/2, r/2), (mid-r/2, r/2)]; continue
        if not (lb > 0):
            nsplit += 1; assert r > arb("1e-9"), ("cannot certify", float(mid.mid()))
            stack += [(mid+r/2, r/2), (mid-r/2, r/2)]; continue
        if ncell % 50000 == 0: print("[progress]",ncell,float(mid.mid()),flush=True)
        if worst is None or lb.lower() < worst[0]: worst = (float(lb.lower()), float(mid.mid()))
        Sl, Sm, Sr = Sc(mid-r), Sc(mid), Sc(mid+r)
        Aq = (Sl+Sr-2*Sm)/(2*r*r); Bq = (Sr-Sl)/(2*r)
        I2 = 2*r*Sm*Sm + 4*Sm*Aq*r**3/3 + 2*Bq*Bq*r**3/3 + 2*Aq*Aq*r**5/5
        e = M3*r**3/(9*s3); Iu = I2.upper(); sq = (2*r*Iu).sqrt() if Iu > 0 else arb(0)
        bint += (I2 + 2*e*sq + 2*r*e*e)/(arb(c) + arb(lb.lower()))
        ssq_low += I2 - 2*e*sq
    T = arb.pi()/2*vnorm2 - arb(ssq_low.lower())
    T = arb(T.upper()) if T.upper() > 0 else arb(0)
    mXB = mlow_tail(XB); assert mXB > 0
    tailb = T/mXB
    beta = (4*(bint + tailb)/arb.pi()).upper()
    margin=(1-C*(arb(beta)+lam*lam/2)).lower()
    coer = (margin*arb(c)).lower()
    out = dict(file=fn, h=f"{hq.p}/{hq.q}", H=f"{Hq.p}/{Hq.q}", N=N, c=f"{c.p}/{c.q}", C=float(C.mid()),
               M_upper=float(Mbern.mid()), Xi=float(Xi.mid()), tail_margin=float(tail.lower()), F2=float(F2.mid()),
               M3=float(M3.mid()), vnorm2=float(vnorm2.mid()), XB=XBi, mlow_XB=float(mXB.lower()),
               parseval_tail_S2=float(T.upper()), beta_tail_part=float((4*tailb/arb.pi()).upper()),
               cells=ncell, splits=nsplit, refinements=nref, REL=REL, prec=ctx.prec,
               worst_cell_lb=worst[0], worst_at=worst[1], beta_upper=float(beta), two_beta=float(2*beta),
               even_coercivity_lower=float(coer), margin_lower=float(margin), lambda_value=float(lam.mid()), verified=bool(margin > 0 and coer > 0), seconds=round(time.time()-t0, 1))
    print(json.dumps(out), flush=True); return out

if __name__ == "__main__":
    res = [verify(fn) for fn in sys.argv[1:]]
    json.dump(res, open(os.environ.get("OUT", "s25_verify_results.json"), "w"), indent=1)

