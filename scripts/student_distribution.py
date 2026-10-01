"""Two-sided Student distribution (regularized incomplete beta)."""
import math

def beta_fraction(a, b, x):
    qab = a + b
    qap = a + 1
    qam = a - 1
    floor = 1e-300
    c = 1.0
    d = 1 - qab * x / qap
    if abs(d) < floor:
        d = floor
    d = 1 / d
    h = d
    for m in range(1, 10001):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d
        c = 1 + aa / c
        if abs(d) < floor:
            d = floor
        if abs(c) < floor:
            c = floor
        d = 1 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d
        c = 1 + aa / c
        if abs(d) < floor:
            d = floor
        if abs(c) < floor:
            c = floor
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < 3e-14:
            return h
    raise ArithmeticError('Incomplete beta did not converge')

def regularized_beta(a, b, x):
    if x == 0:
        return 0.0
    if x == 1:
        return 1.0
    prefactor = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1) / (a + b + 2):
        return prefactor * beta_fraction(a, b, x) / a
    return 1 - prefactor * beta_fraction(b, a, 1 - x) / b

def student_two_sided_p(t, df):
    return regularized_beta(df / 2, 0.5, df / (df + t * t))

def stars(p):
    return '****' if p <= 0.0001 else '***' if p <= 0.001 else '**' if p <= 0.01 else '*' if p <= 0.05 else 'ns'
