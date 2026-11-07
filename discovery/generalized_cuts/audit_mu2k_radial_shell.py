import mpmath as mp

mp.mp.dps = 50

for k in range(1, 7):
    W = lambda r: mp.tanh(r) / mp.cosh(r)**(2*k)
    integral = mp.quad(W, [0, mp.inf])
    exact = mp.mpf(1)/(2*k)
    rstar = mp.atanh(1/mp.sqrt(2*k+1))
    mustar_over_M = mp.mpf('0.5')*mp.sqrt(mp.mpf(2*k)/(2*k+1))
    wstar = 1/mp.sqrt(2*k+1)*(mp.mpf(2*k)/(2*k+1))**k
    print('k=', k)
    print('  integral error =', mp.nstr(abs(integral-exact), 8))
    print('  r*             =', mp.nstr(rstar, 30))
    print('  mu*/M          =', mp.nstr(mustar_over_M, 30))
    print('  W* numeric     =', mp.nstr(W(rstar), 30))
    print('  W* exact       =', mp.nstr(wstar, 30))
    print('  peak error     =', mp.nstr(abs(W(rstar)-wstar), 8))
