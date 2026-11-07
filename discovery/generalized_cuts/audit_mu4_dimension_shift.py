import mpmath as mp

mp.mp.dps = 80


def gamma_ratio(eps):
    return mp.gamma(2 - eps) / mp.gamma(-eps)


def leading_model(eps):
    # Universal UV-pole model I4^(8-2 eps) = Gamma(eps) * Vol(Delta_3)
    return gamma_ratio(eps) * mp.gamma(eps) / 6


if __name__ == '__main__':
    for eps in [mp.mpf('0.1'), mp.mpf('0.03'), mp.mpf('0.01'), mp.mpf('0.003'), mp.mpf('0.001')]:
        ratio = gamma_ratio(eps)
        exact_recursion = -eps * (1 - eps)
        model = leading_model(eps)
        print('eps =', eps)
        print('  gamma-ratio error =', mp.nstr(abs(ratio - exact_recursion), 15))
        print('  UV-residue model  =', mp.nstr(model, 20))
        print('  error to -1/6     =', mp.nstr(abs(model + mp.mpf(1)/6), 15))
