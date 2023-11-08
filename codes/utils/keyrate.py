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

def refresh_l(l,kwargs):
    kwargs['Lgt'] = l
    kwargs['η'] = kwargs['ηd']*10**(-kwargs['ξ']*l/10)

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
        number_of_states = kwargs["number_of_states"]
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
        if mode == 'DP':
            self.PaccjXμ = np.zeros(N)
            self.PaccjXν = np.zeros(N)
            self.PaccjZμ = np.zeros(N)
            self.PaccjZν = np.zeros(N)
            self.Pjμ = np.zeros(N)
            self.Pjν = np.zeros(N)
            for j in range(N):
                self.Pjμ[j] = Pj_β(j, self.μ,kwargs)
                self.Pjν[j] = Pj_β(j, self.ν,kwargs)
                self.PaccjXμ[j] = self.PaccXμ*self.Pjμ[j]
                self.PaccjXν[j] = self.PaccXν*self.Pjν[j]
                self.PaccjZμ[j] = self.PaccZμ*self.Pjμ[j]
                self.PaccjZν[j] = self.PaccZν*self.Pjν[j]
            # Fidelity
            self.FjZμXν = np.zeros(N)
            self.FjZμXμ = np.zeros(N)
            for j in range(N):
                self.FjZμXν[j] = exp(-(self.μ+self.ν)/2)/sqrt(self.Pjμ[j]*self.Pjν[j]) * sum(
                    [(self.μ*self.ν)**((l*N+j)/2)/factorial(l*N+j)**1.5
                        for l in range(number_of_states)])
                self.FjZμXμ[j] = exp(-self.μ)/self.Pjμ[j]*sum(
                    [self.μ**(l*N+j)/factorial(l*N+j)**1.5
                        for l in range(number_of_states)])
            self.Fjμ0 = sqrt(exp(-self.μ)/self.Pjμ[0])
            self.Fjν0 = sqrt(exp(-self.ν)/self.Pjν[0])
            self.FjXμXν = np.zeros(N)
            for j in range(N):
                self.FjXμXν[j] = exp(-(self.μ+self.ν)/2) / \
                    sqrt(self.Pjμ[j]*self.Pjν[j])*sum([(self.μ*self.ν)**((l*N+j)/2)/factorial(l*N+j)
                                                       for l in range(number_of_states)])
            self.FjZμZν = self.FjXμXν
        elif mode == 'CP':
            def P(j, β):
                return exp(-β)*β**j/factorial(j)
            self.Pjμ = np.array([P(j, self.μ) for j in range(cpstates)])
            self.Pjν = np.array([P(j, self.ν) for j in range(cpstates)])
            self.PaccjXμ = self.PaccXμ*self.Pjμ
            self.PaccjXν = self.PaccXν*self.Pjν
            self.PaccjZμ = self.PaccZμ*self.Pjμ
            self.PaccjZν = self.PaccZν*self.Pjν
            self.FjZμXν = np.ones(2)
            self.FjZμXμ = np.ones(2)
            self.FjXμXν = np.ones(cpstates)
            self.Fjμ0 = 1
            self.Fjν0 = 1
            self.FjZμZν = np.ones(cpstates)
        elif mode == 'DP_Cao':
            μ = self.μ
            ν = self.ν

            def P(j, β):
                return exp(-β)*β**j/factorial(j)
            self.Pj = np.zeros(N)
            self.Pjν = np.zeros(N)
            for j in range(N):
                self.Pj[j] = sum([P(l*N+j, self.μ)
                                  for l in range(number_of_states)])
                self.Pjν[j] = sum([P(l*N+j, self.ν)
                                   for l in range(number_of_states)])

            self.Fj = np.zeros(N)

            def Q(m):
                res = 2**(-m/2)*(np.cos(m/4*np.pi)+np.sin(m/4*np.pi))
                return res
            for j in range(N):
                _x = [P(l*N+j, self.μ)*Q(l*N+j)
                      for l in range(number_of_states)]
                self.Fj[j] = sum(_x)/self.Pj[j]
                if self.Fj[j] > 1:
                    self.Fj[j] = 1
            a = sum([(μ*ν)**(l*N/2)/factorial(l*N)
                    for l in range(number_of_states)])
            b = sum([(μ)**(l*N)/factorial(l*N)
                    for l in range(number_of_states)])
            c = sum([(ν)**(l*N)/factorial(l*N)
                    for l in range(number_of_states)])

            self.Fμν = sqrt(a/b/c)
            self.Fμ0 = sqrt(1/b)
            self.Fν0 = sqrt(1/c)

            self.Fμν = np.zeros(N)
            for j in range(N):
                self.Fμν[j] = exp(-(self.μ+self.ν)/2) / \
                    sqrt(self.Pj[j]*self.Pjν[j])*sum([(self.μ*self.ν)**((l*N+j)/2)/factorial(l*N+j)
                                                      for l in range(number_of_states)])
                if self.Fμν[j] > 1:
                    self.Fμν[j] = 1
            '''
            self.Fμ0 = 1
            for i in range(N):
                self.Fj[i] = 1
            self.Fμν = np.ones(N)
            '''

        elif mode == 'CP_Ma':
            self.Qμ = Qeff(self.μ,kwargs)
            self.Qν = Qeff(self.ν,kwargs)
            self.Q0 = Qeff(0,kwargs)
            self.Eμ = Qerr(self.μ,kwargs)/self.Qμ
            self.Eν = Qerr(self.ν,kwargs)/self.Qν
            self.E0 = Qerr(0,kwargs)/self.Q0


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
    if mode == 'DP':
        N = kwargs['N']
    elif mode == 'CP':
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
    mode = kwargs["mode"]
    Ntot = kwargs["Ntot"]
    if mode == 'DP':
        N = kwargs["N"]
    elif mode == 'CP':
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


def Y1_DP_Cao(p: parameters, kwargs):
    N = kwargs["N"]
    nov = 2*N
    A = np.zeros([0, nov])
    b = np.zeros(0)
    for j in range(N):
        a = np.zeros(nov)
        a[j] = 1
        a[N+j] = -1
        A = np.append(A, [a], axis=0)
        b = np.append(b, [sqrt(1-p.Fμν[j]**2)], axis=0)
        A = np.append(A, [-a], axis=0)
        b = np.append(b, [sqrt(1-p.Fμν[j]**2)], axis=0)
    a = np.zeros(nov)
    a[0] = 1
    A = np.append(A, [a], axis=0)
    b = np.append(b, [sqrt(1-p.Fμ0**2)+Qeff(0,kwargs)], axis=0)

    a = np.zeros(nov)
    a[0] = -1
    A = np.append(A, [a], axis=0)
    b = np.append(b, [sqrt(1-p.Fμ0**2)-Qeff(0,kwargs)], axis=0)

    a = np.zeros(nov)
    a[N] = 1
    A = np.append(A, [a], axis=0)
    b = np.append(b, [sqrt(1-p.Fν0**2)+Qeff(0,kwargs)], axis=0)

    a = np.zeros(nov)
    a[N] = -1
    A = np.append(A, [a], axis=0)
    b = np.append(b, [sqrt(1-p.Fν0**2)-Qeff(0,kwargs)], axis=0)

    A1 = np.zeros([0, nov])
    b1 = np.zeros(0)

    a = np.append(p.Pj, np.zeros(N), axis=0)
    A1 = np.append(A1, [a], axis=0)
    b1 = np.append(b1, [Qeff(p.μ,kwargs)], axis=0)

    a = np.append(np.zeros(N), p.Pjν, axis=0)
    A1 = np.append(A1, [a], axis=0)
    b1 = np.append(b1, [Qeff(p.ν,kwargs)], axis=0)

    c = np.zeros(nov)
    c[1] = 1
    bounds = [(0, 1) for i in range(nov)]
    y = [Qeffj(j, p.μ,kwargs) for j in range(N)] + \
        [Qeffj(j, p.ν,kwargs) for j in range(N)]
    '''
    b[0] = sum(y[0:8]*p.Pj)
    b[1] = sum(y[0:8]*p.Pjν)
    '''
    solution = linprog(c, A, b, A1, b1, bounds,
                       'highs-ds')
    Y = solution['x']
    if not Y is None:
        return Y[1]
    else:
        return 0


def e1Y1_DP_Cao(p: parameters, kwargs):
    N = kwargs['N']
    nov = 2*N
    A = np.zeros([0, nov])
    b = np.zeros(0)
    for j in range(N):
        a = np.zeros(nov)
        a[j] = 1
        a[N+j] = -1
        A = np.append(A, [a], axis=0)
        b = np.append(b, [sqrt(1-p.Fμν[j]**2)], axis=0)
        A = np.append(A, [-a], axis=0)
        b = np.append(b, [sqrt(1-p.Fμν[j]**2)], axis=0)

    a = np.zeros(nov)
    a[0] = 1
    A = np.append(A, [a], axis=0)
    b = np.append(b, [sqrt(1-p.Fμ0**2)+Qerr(0,kwargs)], axis=0)

    a = np.zeros(nov)
    a[0] = -1
    A = np.append(A, [a], axis=0)
    b = np.append(b, [sqrt(1-p.Fμ0**2)-Qerr(0,kwargs)], axis=0)

    a = np.zeros(nov)
    a[N] = 1
    A = np.append(A, [a], axis=0)
    b = np.append(b, [sqrt(1-p.Fν0**2)+Qerr(0,kwargs)], axis=0)

    a = np.zeros(nov)
    a[N] = -1
    A = np.append(A, [a], axis=0)
    b = np.append(b, [sqrt(1-p.Fν0**2)-Qerr(0,kwargs)], axis=0)

    A1 = np.zeros([0, nov])
    b1 = np.zeros(0)

    a = np.append(p.Pj, np.zeros(N), axis=0)
    A1 = np.append(A1, [a], axis=0)
    b1 = np.append(b1, [Qerr(p.μ,kwargs)], axis=0)

    a = np.append(np.zeros(N), p.Pjν, axis=0)
    A1 = np.append(A1, [a], axis=0)
    b1 = np.append(b1, [Qerr(p.ν,kwargs)], axis=0)

    c = np.zeros(nov)
    c[1] = 1
    bounds = [(0, 1) for i in range(nov)]
    solution = linprog(-c, A, b, A1, b1, bounds,
                       'highs-ds')
    Y = solution['x']
    if not Y is None:
        return Y[1]
    else:
        return 1


def Y1_CP_Ma(p: parameters):

    μ = p.μ
    ν = p.ν
    Qμ = p.Qμ
    Qν = p.Qν
    Y0 = p.Q0

    Y1 = μ/(μ*ν-ν*ν)*(Qν*exp(ν) - Qμ*exp(μ)*(ν/μ)**2 - (μ**2-ν**2)/μ**2*Y0)
    return Y1


def e1Y1_CP_Ma(p: parameters):
    ν = p.ν
    Qν = p.Qν
    Eν = Qerr(ν)/Qν
    e0Y0 = Qerr(0)
    e1Y1 = (Eν*Qν*exp(ν) - e0Y0)/ν
    return e1Y1


def keyrate(x, l, kwargs):
    refresh_l(l,kwargs)
    if kwargs['mode']=='CPFK':
        K = keyrate_CPFK(kwargs)
        return K.cal(x)
    elif kwargs['mode'] == 'CP_Ma':
        return keyrate_CP_Ma(x,l,kwargs)
    bd = MyBounds()
    isleagal = bd(x_new=np.array(x))
    if not isleagal:
        return -inf
    mode = kwargs["mode"]
    Ntot = kwargs["Ntot"]
    f = kwargs["f"]
    p = parameters(x, kwargs)
    if mode == 'DP_Cao':
        Y1 = Y1_DP_Cao(p, kwargs)
        if Y1 <= 0 or Y1 >= 1:
            return -inf
        e1Y1 = e1Y1_DP_Cao(p,kwargs)
        if e1Y1 >= Y1:
            return -inf
        Δ = (1-p.Fj[1])/(2*Y1)
        if isreal(Δ):
            eb = e1Y1/Y1
            e1 = eb+4*Δ*(1-Δ)*(1-2*eb)+4*(1-2*Δ)*sqrt(absolute(Δ*(1-Δ)*eb*(1-eb)))
        else:
            e1 = 0.5
        R = p.PaccZμ*(p.Pj[1]*Y1*(1-h(e1))-f*Qeff(p.μ,kwargs) *
                      h(Qerr(p.μ,kwargs)/Qeff(p.μ,kwargs)))
        if np.isreal(R):
            return R
        else:
            return -inf
    elif mode == 'CP_Ma':
        #Y1 = self._Y1_CP_Ma()
        # if Y1 == 0:
        #    return -inf
        #e1Y1 = self._e1Y1_CP_Ma()
        #e1 = e1Y1/Y1
        ν = p.ν
        μ = p.μ
        if μ<=ν:
            return -inf
        Y0 = 2*kwargs["pdark"]
        #Qμ = p.Qμ
        #Qν = p.Qν
        H2 = h
        # PRA.72.012326 X.F.Ma eq(36)
        
        # PRA.72.012326 X.F.Ma eq(3)
        η = kwargs["η"]
        Qμ = Y0 + 1-exp(-η*μ)#eq 10
        Qν = Y0 + 1-exp(-η*ν)
        e0 = 1/2
        edet = kwargs["em"]
        Eμ = ((e0*Y0)+edet*(1-exp(-η*μ)))/Qμ #eq 11
        Eν = ((e0*Y0)+edet*(1-exp(-η*ν)))/Qν
        # Q1 = μ**2*exp(-μ)/(μ*ν-ν**2)*(Qν*exp(ν)-Qμ*exp(μ)
        #                               * ν**2/μ**2-(μ**2-ν**2)/μ**2*Y0)#eq35,貌似有问题，不会随着η变化
        Y1Lν0 = μ/(μ*ν-ν**2)*(Qν*exp(ν)-Qμ*exp(μ)*ν**2/μ**2-(μ**2-ν**2)/μ**2*Y0)#eq 34
        Y1 = Y0+1-(1-η)**1#eq6 i=1
        Q1 = Y1*μ*exp(-μ) # eq 8, i = 1
        
        e1 = (Eν*Qν*exp(ν)-e0*Y0)/(Y1Lν0*ν) #eq 37
        # Δ = ν/(μ-ν)*(ν*exp(-ν)*Qμ/(μ*exp(-μ)*Qν)-1)+ν*exp(-ν)*Y0/(μ*Qν)
        q = 0.5
        R = q*(-Qμ*f*H2(Eμ)+Q1*(1-H2(e1)))
        if np.isreal(R):
            return R
        else:
            return -inf
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
