# calculator comes from https://keisan.casio.com/exec/system/14060745333941
from scipy.stats import norm
from numpy import log10
percentiles = {17: 8.493793224109598074445, 18: 8.757290348782315063881,
               19: 9.013271153126674281256, 20: 9.262340089798407573717}
def percentile(n):
    if n > 1e-16:
        return norm.ppf(n)
    elif n < 1e-16:
        f = -int(log10(n))
        if f not in percentiles:
            print('please find the percentile in the given webisite, and put it in the dictionary')
        return percentiles[f]
