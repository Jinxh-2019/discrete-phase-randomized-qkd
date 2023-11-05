# Finite-key analysis for practical implementations of quantum key distribution
from numpy import sqrt, log, exp, log2,Inf
from math import factorial
from codes.utils.entropy import h
from codes.utils.depack import depack_x
import abc
# from codes.utils.keyrate import h, depack_x
# from codes.utils.response_rate import Qeffj, Qerrj
ε_PE = 1e-10
η = 0.1
em = 0.03
pdark = 1e-5
ε_bar = 1e-10
ε_PA = 1e-10
ε_EC = 1e-10
state_length = 64
I = int(1)
II = int(2)
f_EC = 1.1
Ntot = 1e12
# Q = (1-V)/2


class keyrate:
    def __init__(self, kwargs):
        global ε_PE, η, em, pdark, ε_PA, ε_EC, f_EC, Ntot, Q, t
        ε_PE = 1e-10
        η = kwargs['η']
        em = kwargs['em']
        pdark = kwargs['pdark']
        ε_PA = kwargs['εPA']
        ε_EC = kwargs['εcor']
        f_EC = kwargs['f']
        Ntot = kwargs['Ntot']
        # t = kwargs['Lgt']*kwargs['ξ']  # noise
        t=1
        Q = kwargs['Q']

    def cal(self, x):
        global μI, μII, Pμ, Pz_μ, Pν, P0, Pz, Px, N_0, P, f_1_wave, m,R_B,μ
        μI, μII, Pμ, Pz_μ, Pν, P0, Pz, Px = depack_x(x,'CPFK')
        if(μI>μII or μI*exp(-μI)>μII*exp(-μII)): return -Inf
        R_B = [R(0),R(μI),R(μII)]
        μ = [0,μI,μII]
        m = [0, Ntot*(1-Pz_μ)**2, Ntot*0.5**2]
        P = [P0, Pμ, Pν]
        N_0 = P0**2*Ntot
        f_0_wave = R(0)
        f_1_wave = 1/(μII-μI)*(R(I)*μII/p_A(1, I)-R(II)*μI /
                               p_A(1, II))-f_0_wave*(μII+μI)/(μII*μI)
        f_1_wave = min([1,f_1_wave])
        f_1_wave = max([0,f_1_wave])
        if(f_1_wave<0): f_1_wave = 0
        t1 = temp_term(I)
        t2 = temp_term(II)
        if t2>t1: e_X_U_1 = t1  # eq 25
        else: e_X_U_1 = t2
        S_ξ = Y_0_L(I)+Y_1_L(I)*(1-h(e_X_U_1))  # eq 29 S_ξ(A|E,I)
        S_ξ = min([max([0,S_ξ]),1])
        Y0 = 2*pdark
        Qμ = Y0 + 1-exp(-η*μI)
        e0 = 0.5
        Eμ = ((e0*Y0)+em*(1-exp(-η*μI)))/Qμ
        n = Ntot*Pz_μ**2  # P11脚注
        
        Pz_I_2 = (Pz_μ)**2  # p_Z(I)**2
        qI = Pμ  # q(I)
        eZI = e_X(μI)  # e_Z(I)
        # leak_EC = f_EC*h(eZI)+1/n*log2(2/ε_EC)  # eq ?
        leak_EC = f_EC*h(eZI)
        K = qI*R(μI)*Pz_I_2*(S_ξ-Δ(n)-leak_EC)  # eq 29
        return K


def ξ(m, d):
    # see eq3
    part = 2*log(1/ε_PE)+d*log(m+1)
    part /= m
    ans = 1/2*sqrt(part)
    return ans


def f(j):
    Q = ((1-η*em)**(j) - (1-η)**(j)) * \
        (1-pdark)+(1-η)**(j)*(1-pdark)*pdark
    return Q


def Δ(n):
    part1 = log2(2/ε_bar)/n
    part2 = 2/n*log2(1/ε_PA)
    return 7*sqrt(part1)+part2


def p_A(k, γ) -> float:
    μγ = μ[γ]
    return exp(-μγ)*μγ**k/factorial(k)


def R(μ):
    if μ == 0:
        return f(0)
    return 1-(1-2*pdark)*exp(-μ*t*η)


# eq 23
def Y_0_L(γ):
    ans = (p_A(0, γ)*R(0)-ξ(N_0, 2))/R_B[γ]
    if ans<0: return 0
    return ans

# eq 24


def Y_1_L(γ):
    N_γ = Ntot*P[γ]
    ans = (p_A(1, γ)*f_1_wave-ξ(N_γ, 2))/R_B[γ]
    return min([max([0,ans]),1])
# e_X_U(1)


def Y(k, γ):  # eq 27
    term1 = exp(-γ*t*η)
    ans = (1-term1)*Q+term1*pdark
    ans /= R_B[γ]
    return ans


def e_X(μ):
    term1 = exp(-μ*t*η)
    return ((1-term1)*Q+term1*pdark)/R(μ)


def e_X_U(γ):
    if γ == 0:
        return min([1,e_X(μ[γ])+ξ(N_0, 2)])
    return min([1,e_X(μ[γ])+ξ(m[γ], 2)])


def e_X_L(γ):
    if γ == 0:
        return max([e_X(μ[γ])-ξ(N_0, 2),0])
    return max([e_X(μ[γ])-ξ(m[γ], 2),0])


def temp_term(γ): 
    if Y_1_L(γ)==0: return 1
    return (e_X_U(γ)-Y_0_L(γ)*e_X_L(0))/Y_1_L(γ)  # eq 25


# S_ξ(A|E,I)
