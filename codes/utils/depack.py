def depack_x(x, mode):
    if mode not in {}:
        μ = x[0]
        ν = x[1]
        Pμ = x[2]
        Pz_μ = x[3]
        Pν = (1-Pμ)*x[4]
        P0 = (1-Pμ)*(1-x[4])
        Pz = x[2] * x[3] + (1-x[2])*0.5
        Px = 1 - Pz
        return μ, ν, Pμ, Pz_μ, Pν, P0, Pz, Px
    else:
        return None