import mpmath as mp

mp.mp.dps = 80


def shift_factor(r, eps):
    return mp.gamma(r - eps) / mp.gamma(-eps)


def pole_model(n, eps):
    r = n - 2
    return shift_factor(r, eps) * mp.gamma(eps) / mp.factorial(n - 1)


def target(n):
    return -mp.mpf(1) / ((n - 1) * (n - 2))


if __name__ == '__main__':
    for n in range(3, 8):
        print('n =', n, 'target =', mp.nstr(target(n), 20))
        for eps in [mp.mpf('0.03'), mp.mpf('0.01'), mp.mpf('0.003'), mp.mpf('0.001')]:
            value = pole_model(n, eps)
            print('  eps', eps, 'value', mp.nstr(value, 20), 'error', mp.nstr(abs(value-target(n)), 10))
