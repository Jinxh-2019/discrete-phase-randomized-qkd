from random import uniform
def sampling():
    μ = uniform(1e-14, 1-1e-14)
    ν = uniform(1e-14, 1-1e-14)*μ
    Pμ = uniform(1e-14, 1-1e-14)
    fracν0 = uniform(1e-14, 1-1e-14)
    Pz_μ = uniform(1e-14, 1-1e-14)
    x = [μ, ν, Pμ, Pz_μ, fracν0]
    return x