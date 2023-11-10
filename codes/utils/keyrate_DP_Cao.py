import numpy as np
from numpy import sqrt, exp, inf, isreal, absolute
from codes.utils.response_rate import Qeff, Qeffj, Qerr
from math import factorial
# from probabilities import Pj_β
from codes.utils.entropy import h
from scipy.optimize import linprog
from codes.utils.depack import depack_x
def refresh_l(l,kwargs):
    kwargs['Lgt'] = l
    kwargs['η'] = kwargs['ηd']*10**(-kwargs['ξ']*l/10)

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
        self.neffZμ = Ntot*self.PaccZμ*Qeff(self.μ, kwargs)
        self.neffZν = Ntot*self.PaccZν*Qeff(self.ν, kwargs)
        self.nebitXμ = Ntot*self.PaccXμ*Qerr(self.μ, kwargs)
        self.nebitXν = Ntot*self.PaccXν*Qerr(self.ν, kwargs)
        self.neff0 = Ntot*self.P0*Qeff(0, kwargs)
        self.nerr0 = Ntot*self.P0*Qerr(0, kwargs)
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


def keyrate(x, l, kwargs):
    refresh_l(l,kwargs)
    f = kwargs["f"]
    p = parameters(x, kwargs)
    Y1 = Y1_DP_Cao(p, kwargs)
    if Y1 == 0:
        return -inf
    e1Y1 = e1Y1_DP_Cao(p, kwargs)
    Δ = (1-p.Fj[1])/(2*Y1)
    if isreal(Δ):
        eb = e1Y1/Y1
        e1 = eb+4*Δ*(1-Δ)*(1-2*eb)+4*(1-2*Δ)*sqrt(absolute(Δ*(1-Δ)*eb*(1-eb)))
    else:
        e1 = 0.5
    R = p.PaccZμ*(p.Pj[1]*Y1*(1-h(e1))-f*Qeff(p.μ, kwargs) *
                  h(Qerr(p.μ, kwargs)/Qeff(p.μ, kwargs)))
    if np.isreal(R):
        return R
    else:
        return -inf


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
    b = np.append(b, [sqrt(1-p.Fμ0**2)+Qeff(0, kwargs)], axis=0)

    a = np.zeros(nov)
    a[0] = -1
    A = np.append(A, [a], axis=0)
    b = np.append(b, [sqrt(1-p.Fμ0**2)-Qeff(0, kwargs)], axis=0)

    a = np.zeros(nov)
    a[N] = 1
    A = np.append(A, [a], axis=0)
    b = np.append(b, [sqrt(1-p.Fν0**2)+Qeff(0, kwargs)], axis=0)

    a = np.zeros(nov)
    a[N] = -1
    A = np.append(A, [a], axis=0)
    b = np.append(b, [sqrt(1-p.Fν0**2)-Qeff(0, kwargs)], axis=0)

    A1 = np.zeros([0, nov])
    b1 = np.zeros(0)

    a = np.append(p.Pj, np.zeros(N), axis=0)
    A1 = np.append(A1, [a], axis=0)
    b1 = np.append(b1, [Qeff(p.μ, kwargs)], axis=0)

    a = np.append(np.zeros(N), p.Pjν, axis=0)
    A1 = np.append(A1, [a], axis=0)
    b1 = np.append(b1, [Qeff(p.ν, kwargs)], axis=0)

    c = np.zeros(nov)
    c[1] = 1
    bounds = [(0, 1) for i in range(nov)]
    y = [Qeffj(j, p.μ, kwargs) for j in range(N)] + \
        [Qeffj(j, p.ν, kwargs) for j in range(N)]
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
    b = np.append(b, [sqrt(1-p.Fμ0**2)+Qerr(0, kwargs)], axis=0)

    a = np.zeros(nov)
    a[0] = -1
    A = np.append(A, [a], axis=0)
    b = np.append(b, [sqrt(1-p.Fμ0**2)-Qerr(0, kwargs)], axis=0)

    a = np.zeros(nov)
    a[N] = 1
    A = np.append(A, [a], axis=0)
    b = np.append(b, [sqrt(1-p.Fν0**2)+Qerr(0, kwargs)], axis=0)

    a = np.zeros(nov)
    a[N] = -1
    A = np.append(A, [a], axis=0)
    b = np.append(b, [sqrt(1-p.Fν0**2)-Qerr(0, kwargs)], axis=0)

    A1 = np.zeros([0, nov])
    b1 = np.zeros(0)

    a = np.append(p.Pj, np.zeros(N), axis=0)
    A1 = np.append(A1, [a], axis=0)
    b1 = np.append(b1, [Qerr(p.μ, kwargs)], axis=0)

    a = np.append(np.zeros(N), p.Pjν, axis=0)
    A1 = np.append(A1, [a], axis=0)
    b1 = np.append(b1, [Qerr(p.ν, kwargs)], axis=0)

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
