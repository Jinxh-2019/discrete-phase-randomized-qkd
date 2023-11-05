from numpy import log2
def h(x):
    if x >= 0.5:
        return 1
    if x <= 0:
        return 1
    res = -x*log2(x)-(1-x)*log2(1-x)
    return res