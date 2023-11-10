from typing import List, Tuple
from codes.utils.depack import depack_x
from codes.utils.entropy import h
from numpy import sqrt, log, inf, exp, pi, isnan,log2

# system parameters
default = {'pd': 6e-7,
           'p_ap': 0,
           'f': 1.16,
           'em': 5e-3,
           'ε_cor': 1e-15,
           'ε_sec': 1e-8,
           'ηd': 0.1,
           'ξ': 0.2,
           'Ntot': 1e7
           }


def keyrate(x, l, kwargs,ifdefault = False) -> float:
    kwargs_used = {}
    for key in default.keys():
        if not ifdefault and key in kwargs.keys():
            kwargs_used[key] = kwargs[key]
        else:
            kwargs_used[key] = default[key]
    μ, ν, Pμ, Pz_μ, Pν, P0, Pz, Px = depack_x(x, 'DP')
    key_rate = unwrapped_keyrate([Pμ, Pν, μ, ν, Pz], l, kwargs_used)
    return key_rate


def hoeffding(x, y) -> float:
    ans = sqrt(x/2*log(1/y))
    return ans


def unwrapped_keyrate(SeedValue: List, l: int, kwargs=default) -> float:
    n_x = kwargs['Ntot']
    p_dc = kwargs['pd']
    p_ap = kwargs['p_ap']
    f = kwargs['f']
    e_mis = kwargs['em']
    epsilon_cor = kwargs['ε_cor']
    epsilon_sec = kwargs['ε_sec']
    ξ = kwargs['ξ']
    eta_ch = 10**(-ξ*l/10)
    eta_d = kwargs['ηd']
    eta = eta_d*eta_ch
    eps = 21

    p_mu1 = SeedValue[0]
    p_mu2 = SeedValue[1]
    p_mu3 = 1-SeedValue[0]-SeedValue[1]
    mu1 = SeedValue[2]
    mu2 = SeedValue[3]
    mu3 = 0
    q_x = SeedValue[4]
    if max([p_mu1, p_mu2, p_mu3]) > 1 or min([p_mu1, p_mu2, p_mu3])<0 or mu1 >= 1 or mu1 < 0 or mu1 <= mu2+mu3 or mu2 >= mu1 or mu3 >= mu2 or q_x >= 1 or q_x <= 0:
        return -inf
    # 各强度平均探测率（不包括后脉冲贡献）
    D_mu1 = 1-(1-2*p_dc)*exp(-eta*mu1)
    D_mu2 = 1-(1-2*p_dc)*exp(-eta*mu2)
    D_mu3 = 1-(1-2*p_dc)*exp(-eta*mu3)
    # 各强度平均探测率(响应率Q_mui)（包括后脉冲贡献）
    R_mu1 = D_mu1*(1+p_ap)
    R_mu2 = D_mu2*(1+p_ap)
    R_mu3 = D_mu3*(1+p_ap)
    # 各强度比特误码率(E_mu.*Qmui)
    e_mu1 = p_dc+e_mis*(1-exp(-eta*mu1))+p_ap*D_mu1/2
    e_mu2 = p_dc+e_mis*(1-exp(-eta*mu2))+p_ap*D_mu2/2
    e_mu3 = p_dc+e_mis*(1-exp(-eta*mu3))+p_ap*D_mu3/2
    # Alice发送总脉冲数
    N = n_x/(q_x*q_x)/(p_mu1*R_mu1+p_mu2*R_mu2+p_mu3*R_mu3)
    # x基各强度样本数
    n_x_mu1 = N*(q_x*q_x)*p_mu1*R_mu1
    n_x_mu2 = N*(q_x*q_x)*p_mu2*R_mu2
    n_x_mu3 = N*(q_x*q_x)*p_mu3*R_mu3
    # z基各强度样本数
    n_z = (1-q_x)*(1-q_x)*n_x/q_x/q_x  # Alice发送Z基的个数
    n_z_mu1 = N*((1-q_x)*(1-q_x))*p_mu1*R_mu1
    n_z_mu2 = N*((1-q_x)*(1-q_x))*p_mu2*R_mu2
    n_z_mu3 = N*((1-q_x)*(1-q_x))*p_mu3*R_mu3

    tau_0 = p_mu1*exp(-mu1)+p_mu2*exp(-mu2)+p_mu3*exp(-mu3)

    # x基真空事件数
    n_x_mu3_l = max(
        [0, exp(mu3)*(n_x_mu3-hoeffding(n_x, epsilon_sec/eps))/p_mu3])
    n_x_mu2_u = exp(mu2)*(n_x_mu2+hoeffding(n_x, epsilon_sec/eps))/p_mu2
    s_x_0_l = max([0, tau_0*(mu2*n_x_mu3_l-mu3*n_x_mu2_u)/(mu2-mu3)])

    # z基真空事件数
    n_z_mu3_l = max(
        [0, exp(mu3)*(n_z_mu3-hoeffding(n_z, epsilon_sec/eps))/p_mu3])
    n_z_mu2_u = exp(mu2)*(n_z_mu2+hoeffding(n_z, epsilon_sec/eps))/p_mu2
    s_z_0_l = max([0, tau_0*(mu2*n_z_mu3_l-mu3*n_z_mu2_u)/(mu2-mu3)])

    tau1 = p_mu1*exp(-mu1)*mu1+p_mu2*exp(-mu2)*mu2+p_mu3*exp(-mu3)*mu3

    # x基单光子事件数
    n_x_mu2_l = max(
        [0, exp(mu2)*(n_x_mu2-hoeffding(n_x, epsilon_sec/eps))/p_mu2])
    n_x_mu3_u = exp(mu3)*(n_x_mu3+hoeffding(n_x, epsilon_sec/eps))/p_mu3
    n_x_mu1_u = exp(mu1)*(n_x_mu1+hoeffding(n_x, epsilon_sec/eps))/p_mu1
    s_x_1_l = max([0, tau1*mu1/(mu1*(mu2-mu3)-mu2*mu2+mu3*mu3)*(n_x_mu2_l -
                  n_x_mu3_u-((mu2 ** 2-mu3 ** 2)/mu1 ** 2)*(n_x_mu1_u-s_x_0_l/tau_0))])

    #z基单光子事件数
    n_z_mu2_l=max([0,exp(mu2)*(n_z_mu2-hoeffding(n_z,epsilon_sec/eps))/p_mu2])
    n_z_mu3_u=exp(mu3)*(n_z_mu3+hoeffding(n_z,epsilon_sec/eps))/p_mu3
    n_z_mu1_u=exp(mu1)*(n_z_mu1+hoeffding(n_z,epsilon_sec/eps))/p_mu1
    s_z_1_l=max([0,tau1*mu1/(mu1*(mu2-mu3)-mu2*mu2+mu3*mu3)*(n_z_mu2_l-n_z_mu3_u-((mu2**2-mu3**2)/mu1**2)*(n_z_mu1_u-s_z_0_l/tau_0))])
    
    #z基单光子比特误码数
    m_z=N*(1-q_x)*(1-q_x)*(p_mu1*e_mu1+p_mu2*e_mu2+p_mu3*e_mu3) #Z基总误码数
    m_z_mu2=N*(1-q_x)*(1-q_x)*p_mu2*e_mu2
    m_z_mu3=N*(1-q_x)*(1-q_x)*p_mu3*e_mu3
    m_z_mu2_u=exp(mu2)/p_mu2*(m_z_mu2+hoeffding(m_z,epsilon_sec/eps))
    m_z_mu3_l=max(0,exp(mu3)/p_mu3*(m_z_mu3-hoeffding(m_z,epsilon_sec/eps)))
    v_z_1_u=tau1*(m_z_mu2_u-m_z_mu3_l)/(mu2-mu3) #Z基单光子误码数
    
    #x基单光子事件相位误码率
    gamma=RanSam(epsilon_sec/eps,v_z_1_u/s_z_1_l,s_z_1_l,s_x_1_l)
    phi_x_u=v_z_1_u/s_z_1_l+gamma

    #x基总误码率
    m_x=N*q_x**2*(p_mu1*e_mu1+p_mu2*e_mu2+p_mu3*e_mu3) #X基总误码数
    e_obs=m_x/n_x

    ll=s_x_0_l+s_x_1_l-s_x_1_l*h(phi_x_u)-f*h(e_obs)*n_x-6*log2(eps/epsilon_sec)-log2(2/epsilon_cor)
    R=ll/N
    return R
def RanSam(a,b,c,d):
    gamma = sqrt((c+d)*(1-b)*b/(c*d)*log((c+d)/(2*pi)/(a*a*(c*d*(1-b)*b))))
    if gamma<0 or isnan(gamma):
        return 0
    return gamma