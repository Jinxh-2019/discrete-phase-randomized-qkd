from numpy import exp
from math import factorial
from codes.utils.probabilities import Pj_β


def Qcorr(μ, kwargs):  # 计算比特误码用到的函数, μ is float
    # self.update()
    #        μ = μ.β
    # pdark = self.pdark
    # η = self.η
    # em = self.em
    pdark = kwargs['pdark']
    η = kwargs['η']
    em = kwargs['em']
    res = (1-(1-pdark)*exp(-η*(1-em)*μ))*exp(-η*em*μ)*(1-pdark)
    return res


def Qerr(μ, kwargs):  # 计算比特误码用到的函数,μ is float
    #        μ = μ.β
    # self.update()
    pdark = kwargs['pdark']
    η = kwargs['η']
    em = kwargs['em']
    res = (1-pdark)*exp(-η*(1-em)*μ)*(1-(1-pdark)*exp(-η*em*μ))
    return res


def Qeff(μ, kwargs):
    return Qcorr(μ, kwargs)+Qerr(μ, kwargs)

# def Qcorri(i, μ, kwargs):
#     η = kwargs['η']
#     em = kwargs['em']
#     pdark = kwargs['pdark']
#     Q = ((1-η*em)**(i) - (1-η)**(i)) * \
#                 (1-pdark)+(1-η)**(i)*(1-pdark)*pdark
def Qcorrj(j, μ, kwargs):
    pdark = kwargs['pdark']
    η = kwargs['η']
    em = kwargs['em']
    mode = kwargs['mode']
    N = kwargs['N']
    def f(l):
        Q = ((1-η*em)**(l*N+j) - (1-η)**(l*N+j)) * \
            (1-pdark)+(1-η)**(l*N+j)*(1-pdark)*pdark
        return Q
    return f(j)


def Qerrj(j, μ, kwargs):  # 计算比特误码用到的函数,μ is float
    #        μ = μ.β
    # self.update()
    pdark = kwargs['pdark']
    η = kwargs['η']
    em = kwargs['em']
    mode = kwargs['mode']
    N = kwargs['N']
    def f(l):
        Q = ((1-η*(1-em))**(l*N+j)-(1-η)**(l*N+j)) * \
            (1-pdark)+(1-η)**(l*N+j)*(1-pdark)*pdark
        return Q
    return f(j)

def Qeffj(j, μ, kwargs):
        return Qcorrj(j, μ, kwargs)+Qerrj(j, μ, kwargs)