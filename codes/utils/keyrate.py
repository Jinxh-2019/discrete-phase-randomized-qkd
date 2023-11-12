from codes.utils.response_rate import Qerr, Qeff, Qcorr, Qeffj, Qerrj
from codes.utils.probabilities import Pj_β
from codes.utils.keyrate_CPFK import keyrate as keyrate_CPFK
from codes.utils.entropy import h
from codes.utils.depack import depack_x
import numpy as np
from codes.utils.lemma1 import Δ
from numpy import exp, sqrt, log2, inf,absolute,isreal
from math import factorial,isnan
from scipy.optimize import linprog
from codes.utils.keyrate_CP_Ma import keyrate as keyrate_CP_Ma
from codes.utils.keyrate_DPFK import keyrate as keyrate_DPFK
from codes.utils.keyrate_DP_Cao import keyrate as keyrate_DP_Cao
from codes.utils.keyrate_CPFK_mine import keyrate as keyrate_CPFK_mine
from codes.utils.keyrate_decoy_Wang import keyrate as keyrate_decoy

def keyrate(x, l, kwargs):
    mode = kwargs['mode']
    if mode == 'DP_Cao':
        return keyrate_DP_Cao(x, l, kwargs)
    elif mode == 'CP_Ma':
        return keyrate_CP_Ma(x, l, kwargs)
    elif mode == 'DPFK':
        return keyrate_DPFK(x, l, kwargs)
    elif mode == 'CPFK':
        return keyrate_CPFK(x, l, kwargs)
    elif mode == 'CPFK_mine':
        return keyrate_CPFK_mine(x, l, kwargs)
