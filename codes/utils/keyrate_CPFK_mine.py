import numpy as np
from numpy import sqrt, exp, log2, inf, isnan
from math import factorial
from codes.utils.entropy import h
from codes.utils.depack import depack_x
from codes.utils.response_rate_continuouse_phase_randomized import Qeff, Qeffj, Qerr
from codes.utils.probabilities import Pj_β
from codes.utils.lemma1 import Δ
from scipy.optimize import linprog
from decimal import Decimal, getcontext, localcontext
def refresh_l(l,kwargs):
    kwargs['Lgt'] = l
    kwargs['η'] = kwargs['ηd']*10**(-kwargs['ξ']*l/10)

def keyrate(x, l, kwargs):
    refresh_l(l,kwargs)
    bd = MyBounds()
    isleagal = bd(x_new=np.array(x))
    if not isleagal:
        return -inf
    # mode = kwargs["mode"]
    Ntot = kwargs["Ntot"]
    f = kwargs["f"]
    p = parameters(x, kwargs)
    f = kwargs['f']
    εcor = kwargs["εcor"]
    εPA = kwargs["εPA"]
    _ebit = ebit(x[0], kwargs)
    _eph = eph(p, kwargs)
    _n1 = p.n1
    nbit = p.neffZμ
    keylength = _n1*(1-h(_eph)) - nbit*f*h(_ebit) - \
        (log2(2/εcor)+log2(1/εPA ** 2/4))
    keyrate = keylength/Ntot
    # print(keyrate,x)
    return keyrate

class MyBounds:
    def __init__(self, xmax=[1-1e-14, 1-1e-14, 1-1e-14, 1-1e-14, 1-1e-14], xmin=[1e-14, 1e-14, 1e-14, 1e-14, 1e-14]):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)

    def __call__(self, x_new):
        x = x_new
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin
class parameters:
    def __init__(self, x, kwargs):
        N = kwargs["N"]
        Ntot = kwargs["Ntot"]
        mode = kwargs["mode"]
        # number_of_states = kwargs["number_of_states"]
        cpstates = kwargs["cpstates"]
        self.is_n1 = False
        self.x = x
        self.μ, self.ν, self.Pμ, self.Pz_μ, self.Pν, self.P0, self.Pz, self.Px = depack_x(
            x, mode)
        self.PaccZν = self.Pν*0.5*self.Pz
        self.PaccXν = self.Pν*0.5*self.Px
        self.PaccZμ = self.Pμ*self.Pz_μ*self.Pz
        self.PaccXμ = self.Pμ*(1-self.Pz_μ)*self.Px
        # Counts
        self.neffZμ = Ntot*self.PaccZμ*Qeff(self.μ,kwargs)
        self.neffZν = Ntot*self.PaccZν*Qeff(self.ν,kwargs)
        self.nebitXμ = Ntot*self.PaccXμ*Qerr(self.μ,kwargs)
        self.nebitXν = Ntot*self.PaccXν*Qerr(self.ν,kwargs)
        self.neff0 = Ntot*self.P0*Qeff(0,kwargs)
        self.nerr0 = Ntot*self.P0*Qerr(0,kwargs)
        def P(j, β):
            return exp(-β)*β**j/factorial(j)
        self.Pjμ = np.array([P(j, self.μ) for j in range(cpstates)])
        self.Pjν = np.array([P(j, self.ν) for j in range(cpstates)])
        self.PaccjXμ = self.PaccXμ*self.Pjμ
        self.PaccjXν = self.PaccXν*self.Pjν
        self.PaccjZμ = self.PaccZμ*self.Pjμ
        self.PaccjZν = self.PaccZν*self.Pjν
        with localcontext() as ctx:
            ctx.prec = 100
            self.FjZμXν = [Decimal('1')]*2
            self.FjZμXμ = [Decimal('1')]*2
            self.FjXμXν = [Decimal('1')]*cpstates
            self.Fjμ0 = Decimal('1')
            self.Fjν0 = Decimal('1')
            self.FjZμZν = [Decimal('1')]*cpstates

def ebit(μ, kwargs):
    Qe = Qerr(μ, kwargs)
    Qef = Qeff(μ, kwargs)
    return Qe/Qef

def constraint_type1(A, b, a1, a2, P1, P2, F, n1, n2, nov, kwargs):
    a = np.zeros(nov)
    if P1 > P2:
        a[a1] = P2/P1
        a[a2] = -1
        A = np.append(A, [a], axis=0)
        b = np.append(b, [Δ(P1, P2, F, n1, n2, kwargs)], axis=0)
        if isnan(b[-1]):
            print()
    else:
        a[a1] = 1
        a[a2] = -P1/P2
        A = np.append(A, [a], axis=0)
        b = np.append(b, [Δ(P2, P1, F, n2, n1, kwargs)], axis=0)
        if isnan(b[-1]):
            print()
    a = np.zeros(nov)
    if P1 > P2:
        a[a1] = -P2/P1
        a[a2] = 1
        A = np.append(A, [a], axis=0)
        b = np.append(b, [b[-1]], axis=0)
    else:
        a[a1] = -1
        a[a2] = P1/P2
        A = np.append(A, [a], axis=0)
        b = np.append(b, [b[-1]], axis=0)
    return A, b


def constraint_type2(A, b, a1, P1, P2, F, n1, n2, nov, kwargs):
    a = np.zeros(nov)
    if P1 > P2:
        a[a1] = P2/P1
        delta = Δ(P1, P2, F, n1, n2, kwargs)
        if isnan(delta):
            print()
        A = np.append(A, [a], axis=0)
        b = np.append(b, [delta+n2], axis=0)
    else:
        a[a1] = 1
        delta = Δ(P2, P1, F, n2, n1, kwargs)
        if isnan(delta):
            print()
        A = np.append(A, [a], axis=0)
        b = np.append(b, [delta+P1/P2*n2], axis=0)
    a = np.zeros(nov)
    if P1 > P2:
        a[a1] = -P2/P1
        A = np.append(A, [a], axis=0)
        b = np.append(b, [delta-n2], axis=0)
    else:
        a[a1] = -1
        A = np.append(A, [a], axis=0)
        b = np.append(b, [delta-P1/P2*n2], axis=0)
    return A, b

def n1(p: parameters, kwargs):
    if p.is_n1:
        return p.n1
    Ntot = kwargs["Ntot"]
    mode = kwargs['mode']
    N = kwargs['cpstates']
    nov = N*2
    A = np.zeros([0, nov])
    b = np.zeros(0)
    for j in range(N):
        P1 = p.PaccjZμ[j]
        P2 = p.PaccjZν[j]
        F = p.FjZμZν[j]
        n1 = p.neffZμ
        n2 = p.neffZν
        A, b = constraint_type1(A, b, j, N+j, P1, P2, F, n1, n2, nov, kwargs)

    P1 = p.PaccjZμ[0]
    P2 = p.P0
    F = p.Fjμ0
    n1 = p.neffZμ
    n2 = p.neff0
    A, b = constraint_type2(A, b, 0, P1, P2, F, n1, n2,nov, kwargs)
    c = np.zeros(nov)

    P1 = p.PaccjZν[0]
    P2 = p.P0
    F = p.Fjν0
    n1 = p.neffZμ
    n2 = p.neff0

    A, b = constraint_type2(A, b, N, P1, P2, F, n1, n2,nov, kwargs)
    c = np.zeros(nov)
    c[1] = 1
    bounds = [(0, p.PaccjZμ[i]*Ntot) for i in range(N)] + \
        [(0, p.PaccjZν[i]*Ntot) for i in range(N)]
    A1 = np.zeros([0, nov])
    b1 = np.zeros(0)

    a = [1 for i in range(N)]+[0 for i in range(N)]
    A1 = np.append(A1, [a], axis=0)
    b1 = np.append(b1, [p.neffZμ], axis=0)

    a = [0 for i in range(N)]+[1 for i in range(N)]
    A1 = np.append(A1, [a], axis=0)
    b1 = np.append(b1, [p.neffZν], axis=0)

    solution = linprog(c, A, b, A1, b1, bounds,
                       'highs-ds')
    if 'x' in solution and solution['x'] is not None:
        n1 = solution['x'][1]
    else:
        n1 = 0
    p.n1 = n1
    return n1

def neph1(p: parameters, kwargs):
    Ntot = kwargs["Ntot"]
    N = kwargs["cpstates"]
    nov = N*2+1
    A = np.zeros([0, nov])
    b = np.zeros(0)
    p.nbit = p.neffZμ
    # 6
    for j in range(N):
        P1 = p.PaccjXμ[j]
        P2 = p.PaccjXν[j]
        F = p.FjXμXν[j]
        n1 = p.nebitXμ
        n2 = p.nebitXν
        A, b = constraint_type1(A, b, j, N+j, P1, P2, F, n1, n2, nov, kwargs)
    # 5
    P1 = p.PaccjXμ[0]
    P2 = p.P0
    F = p.Fjμ0
    n1 = p.nebitXμ
    n2 = p.nerr0
    A, b = constraint_type2(A, b, 0, P1, P2, F, n1, n2, nov, kwargs)
    # 5
    P1 = p.PaccjXν[0]
    P2 = p.P0
    F = p.Fjν0
    n1 = p.nebitXν
    n2 = p.nerr0
    A, b = constraint_type2(A, b, N, P1, P2, F, n1, n2, nov, kwargs)
    # 4
    P1 = p.PaccjZμ[1]
    P2 = p.PaccjXμ[1]
    F = p.FjZμXμ[1]
    n1 = p.neffZμ
    n2 = p.nebitXμ
    A, b = constraint_type1(A, b, 2*N, 1, P1, P2, F, n1, n2, nov, kwargs)
    # 3
    P1 = p.PaccjZμ[1]
    P2 = p.PaccjXν[1]
    F = p.FjZμXν[1]
    n1 = p.neffZμ
    n2 = p.nebitXν
    A, b = constraint_type1(A, b, 2*N, 1+N, P1, P2, F, n1, n2, nov, kwargs)

    c = np.zeros(nov)
    c[-1] = 1
    bounds = [(0, p.PaccjXμ[i]*Ntot) for i in range(N)] + \
        [(0, p.PaccjXν[i]*Ntot) for i in range(N)]+[(0, None)]
    A1 = np.zeros([0, nov])
    b1 = np.zeros(0)

    a = [1 for i in range(N)]+[0 for i in range(N)]+[0]
    A1 = np.append(A1, [a], axis=0)
    b1 = np.append(b1, [p.nebitXμ], axis=0)

    a = [0 for i in range(N)]+[1 for i in range(N)]+[0]
    A1 = np.append(A1, [a], axis=0)
    b1 = np.append(b1, [p.nebitXν], axis=0)

    solution = linprog(-c, A, b, A1, b1, bounds,
                       'highs-ds')
    if 'x' in solution and solution['x'] is not None:
        n1 = solution['x'][-1]
    else:
        n1 = 0
    return n1

def eph(p: parameters, kwargs):
    _n1 = n1(p, kwargs)
    if _n1 == 0:
        return 1
    return neph1(p, kwargs)/_n1
