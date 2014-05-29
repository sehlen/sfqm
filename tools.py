from sage.all import RealField, floor, pari, cached_function
from psage.modules.finite_quadratic_module import FiniteQuadraticModule

def BB(x):
    RF = RealField(100)
    x = RF(x)
    onehalf = RF(1) / 2
    return x - onehalf * (floor(x) - floor(-x))


@cached_function
def h(d, prec=53):
    RF = RealField(prec)
    kappa = lambda d: RF(1) / 3 if d == 3 else (RF(1) / 2 if d == 4 else 1)
    return RF(pari("qfbclassno(%s,1)" % (-d))) * kappa(d)
