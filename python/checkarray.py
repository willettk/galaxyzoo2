import numpy as N

def checkarray(x, double=False):
    if double:
        t = N.float64
    else:
        t = N.float32
    xa = N.asarray(x, t)
    if len(xa.shape) < 1:
        xa = N.asarray([x], t)
    return xa
