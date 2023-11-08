from numpy import sqrt, log
from codes.utils.calculated_percentiles import percentile
from scipy.stats import norm


def δ(x, y, z):
    if x < 0:
        return 0
    else:
        res = sqrt(3*x*y*log(1/z))
        return res


def δ_1(n, p, e):
    if n < 0:
        return 0
    μ = n*p
    if μ < 5:
        return δ(n, p, e)
    σ = sqrt(n*p*(1-p))
    ppf = percentile(e)
    return ppf*σ


def N1(P, kwargs):
    Ntot = kwargs['Ntot']
    ε1 = kwargs['ε1']
    d = δ(Ntot, P, ε1)
    N1 = Ntot*P + d
    return N1


def N2(P, kwargs):
    Ntot = kwargs['Ntot']
    ε1 = kwargs['ε1']
    d = δ(Ntot, P, ε1)
    N2 = Ntot*P - d
    return max(N2, 0)


def Δ(P1, P2, F, n1, n2, kwargs):
    if F > 1:
        F = 1
    lemma = kwargs['lemma']
    if lemma == 'lemmaA1':
        return Δ_lemmaA1(P1, P2, F, n1, n2, kwargs)
    elif lemma == 'new_lemma':
        return Δ_new(P1, P2, F, kwargs)
    elif lemma == 'lemmaA1_without_the_first_term':
        return Δ_lemmaA1_1(P1, P2, F, n1, n2, kwargs)
    elif lemma == 'lemmaA1_without_the_second_term':
        return Δ_lemmaA1_2(P1, P2, F, n1, n2, kwargs)
    elif lemma == 'lemmaA1_without_the_third_term':
        return Δ_lemmaA1_3(P1, P2, F, n1, n2, kwargs)
    elif lemma == 'lemma_shan':
        return Δ1emma_shan(P1, P2, F, n1, n2, kwargs)
    elif lemma == 'lemma_shan2':
        return Δ1emma_shan2(P1, P2, F, n1, n2, kwargs)
    elif lemma == 'F_test':
        return Δ_lemmaA1(P1, P2, 1, n1, n2, kwargs)


def Δ_new(P1, P2, F, kwargs):
    N = (P1+P2)*kwargs['Ntot']
    ε = kwargs['ε0']
    εμ = kwargs['ε1']
    rN = sqrt(N)
    Δ = rN/2*norm.cdf(1-ε)+sqrt(1-F*F)*N+rN/2*norm.cdf(1-εμ)
    return Δ


def Δ_lemmaA1(P1, P2, F, n1, n2, kwargs):
    _N1 = N1(P2, kwargs)
    _N2 = N2(P2, kwargs)
    ε0 = kwargs['ε0_A1']
    ε1 = kwargs['ε1_A1']
    ε2 = kwargs['ε2_A1']
    if F > 1:
        F = 1
    Δ = _N1*sqrt(1-F*F) + 2*δ(_N1, (1+sqrt(1-F*F))/2, ε0) - \
        2*δ(_N2-n1-n2, 1/2, ε1) + δ(n1, P2/P1, ε2)
    #Δ = N1*sqrt(1-F*F)
    return Δ


def Δ_lemmaA1_1(P1, P2, F, n1, n2, kwargs):
    _N1 = N1(P2, kwargs)
    _N2 = N2(P2, kwargs)
    ε0 = kwargs['ε0_A1']
    ε1 = kwargs['ε1_A1']
    ε2 = kwargs['ε2_A1']
    if F > 1:
        F = 1
    Δ = _N1*sqrt(1-F*F) - 2*δ(_N2-n1-n2, 1/2, ε1) + δ(n1, P2/P1, ε2)
    #Δ = N1*sqrt(1-F*F)
    return Δ


def Δ_lemmaA1_2(P1, P2, F, n1, n2, kwargs):
    _N1 = N1(P2, kwargs)
    _N2 = N2(P2, kwargs)
    ε0 = kwargs['ε0_A1']
    ε1 = kwargs['ε1_A1']
    ε2 = kwargs['ε2_A1']
    if F > 1:
        F = 1
    Δ = _N1*sqrt(1-F*F) + 2*δ(_N1, (1+sqrt(1-F*F))/2, ε0) + δ(n1, P2/P1, ε2)
    return Δ


def Δ_lemmaA1_3(P1, P2, F, n1, n2, kwargs):
    _N1 = N1(P2, kwargs)
    _N2 = N2(P2, kwargs)
    ε0 = kwargs['ε0_A1']
    ε1 = kwargs['ε1_A1']
    ε2 = kwargs['ε2_A1']
    if F > 1:
        F = 1
    Δ = _N1*sqrt(1-F*F) + 2*δ(_N1, (1+sqrt(1-F*F))/2, ε0) - \
        2*δ(_N2-n1-n2, 1/2, ε1)
    return Δ
# def Δ1emma_shan1(P1, P2, F, n1, n2, kwargs):
#     _N1 = N1(P2, kwargs)
#     _N2 = N2(P2, kwargs)
#     ε0 = kwargs['ε0']
#     ε1 = kwargs['ε1']
#     ε2 = kwargs['ε2']
#     if F > 1:
#         F = 1
#     Δ = _N1*sqrt(1-F*F)/2 + 2*δ(_N1, (1+sqrt(1-F*F))/2, ε0) - \
#         2*δ(_N2-n1-n2, 1/2, ε1) + δ(n1, P2/P1, ε2)
#     #Δ = N1*sqrt(1-F*F)
#     return Δ


def Δ1emma_shan(P1, P2, F, n1, n2, kwargs):
    Ntot = kwargs['Ntot']
    ε0 = kwargs['ε0']
    ε1 = kwargs['ε1']
    ε2 = kwargs['ε2']
    _N1 = N1(P2, kwargs)
    _N2 = N2(P2, kwargs)
    if F > 1:
        F = 1
    Δ = Ntot*P2*sqrt(1-F*F)/2 + Ntot*P2-_N2 + δ(P2/P1,n1,ε2)
    
    #Δ = N1*sqrt(1-F*F)
    return Δ


def Δ1emma_shan2(P1, P2, F, n1, n2, kwargs):
    _N1 = N1(P2, kwargs)
    _N2 = N2(P2, kwargs)
    ε0 = kwargs['ε0']
    ε1 = kwargs['ε1']
    ε2 = kwargs['ε2']
    if F > 1:
        F = 1
    Δ = _N1*sqrt(1-F*F)/2 + 2*δ(_N1, (1+sqrt(1-F*F))/2, ε0) - \
        2*δ(_N2-n1-n2, 1/2, ε1) + δ(n1, P2/P1, ε2)
    #Δ = N1*sqrt(1-F*F)
    return Δ

# def Δ1_4(P1, P2, F, n1, n2, kwargs):
#     Ntot = kwargs['Ntot']
#     _N1 = N1(P2, kwargs)
#     _N2 = N2(P2, kwargs)
#     ε0 = kwargs['ε0']
#     ε1 = kwargs['ε1']
#     ε2 = kwargs['ε2']
#     if F > 1:
#         F = 1
#     Δ = Ntot*P2*sqrt(1-F*F) + δ_1(Ntot, P2*(1+sqrt(1-F*F))/2, ε0) - \
#         δ_1(n1, P2/P1, ε1) + δ_1(Ntot, P2, ε2)
#     #Δ = N1*sqrt(1-F*F)
#     return Δ

# def Δ1_1_1(P1, P2, F, n1, n2, kwargs):
#     _N1 = N1(P2, kwargs)
#     _N2 = N2(P2, kwargs)
#     ε0 = kwargs['ε0']
#     ε1 = kwargs['ε1']
#     ε2 = kwargs['ε2']
#     if F > 1:
#         F = 1
#     Δ = _N1*sqrt(1-F*F)/2 + 2*δ(_N1, (1+sqrt(1-F*F))/2, ε0) + \
#         2*δ(_N2-n1-n2, 1/2, ε1) + δ(n1, P2/P1, ε2)
#     #Δ = N1*sqrt(1-F*F)
#     return Δ

# def Δ1_0_2(P1, P2, F, n1, n2, kwargs):
#     _N1 = N1(P2, kwargs)
#     _N2 = N2(P2, kwargs)
#     ε0 = kwargs['ε0']
#     ε1 = kwargs['ε1']
#     ε2 = kwargs['ε2']
#     if F > 1:
#         F = 1
#     Δ = _N1*sqrt(1-F*F) + 2*δ(_N1, (1+sqrt(1-F*F))/2, ε0) - \
#         2*δ(_N2-n1-n2, 1/2, ε1)*0 + δ(n1, P2/P1, ε2)
#     #Δ = N1*sqrt(1-F*F)
#     return Δ

# def Δ2(P1, P2, F, n1, n2, kwargs):
#     _N1 = N1(P2, kwargs)
#     _N2 = N2(P2, kwargs)
#     ε0 = kwargs['ε0']
#     ε1 = kwargs['ε1']
#     ε2 = kwargs['ε2']
#     if F > 1:
#         F = 1
#     Δ = _N1*sqrt(1-F*F) + 2*δ(_N1, (1+sqrt(1-F*F))/2, ε0) - \
#         2*δ(_N2-n1-n2, 1/2, ε1)*0 + δ(n1, P2/P1, ε2)
#     #Δ = N1*sqrt(1-F*F)
#     return Δ
