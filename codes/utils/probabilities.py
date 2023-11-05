from numpy import exp
from math import factorial
def Pj_β(j, β,kwargs):
    number_of_states = kwargs['number_of_states']
    N = kwargs['N']
    Pj_β = exp(-β)*sum([β**(l*N+j)/factorial(l*N+j)
                        for l in range(number_of_states)])
    return Pj_β