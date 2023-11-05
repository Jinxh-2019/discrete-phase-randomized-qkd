
Y1 = Y1_DP_Cao(p, kwargs)
if Y1 == 0:
    return -inf
e1Y1 = e1Y1_DP_Cao(p,kwargs)
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