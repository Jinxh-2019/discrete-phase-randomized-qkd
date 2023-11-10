from numpy import exp,inf,isreal
from codes.utils.entropy import h
from codes.utils.depack import depack_x

def refresh_l(l,kwargs):
    kwargs['Lgt'] = l
    kwargs['η'] = kwargs['ηd']*10**(-kwargs['ξ']*l/10)
    return kwargs['η']

def Q(μ,kwargs):
    Y0 = 2*kwargs['pdark']
    ans = Y0+1-exp(-kwargs['η']*μ)
    return ans
def keyrate(x,l,kwargs):
    η = refresh_l(l,kwargs)
    f = kwargs['f']
    #Y1 = self._Y1_CP_Ma()
    # if Y1 == 0:
    #    return -inf
    #e1Y1 = self._e1Y1_CP_Ma()
    #e1 = e1Y1/Y1
    μ, ν, Pμ, Pz_μ, Pν, P0, Pz, Px=depack_x(x,'CP_Ma')
    Y0 = 2*kwargs['pdark']
    
    if μ<=ν:
        return -inf
    Y0 = 2*kwargs["pdark"]
    #Qμ = p.Qμ
    #Qν = p.Qν
    H2 = h
    # PRA.72.012326 X.F.Ma eq(36)

    # PRA.72.012326 X.F.Ma eq(3)
    Qμ = Y0+1-exp(-η*μ)
    Qν = Y0+1-exp(-η*μ)
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
    q = 1
    R = q*(-Qμ*f*H2(Eμ)+Q1*(1-H2(e1)))
    if isreal(R):
        return R
    else:
        return -inf

