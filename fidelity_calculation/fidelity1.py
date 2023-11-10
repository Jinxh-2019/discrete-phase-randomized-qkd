import matplotlib.pyplot as plt
from decimal import Decimal, getcontext, localcontext
from math import factorial as math_factorial

# Set the precision (for example, 50 decimal places)
getcontext().prec = 50


def factorial(n):
    # We use the standard math factorial but convert the result to a Decimal
    return Decimal(math_factorial(n))


# Assuming μ and ν are already Decimal objects and Pjμ, Pjν are lists of Decimal objects
# and number_of_states is an integer
Y = []
for N in range(2,32,2):

    number_of_states = 2
    FjZμXν = [0] * N
    FjXμXν = [0] * N
    with localcontext() as ctx:  # Use local context to avoid affecting global precision
        ctx.prec = 100  # Set precision for this block
        μ = Decimal('0.5')
        ν = Decimal('0.005')
        Pjμ = [(-μ).exp()*sum([μ**(l*N+j)/factorial(l*N+j) for l in range(number_of_states)]
                            )for j in range(N)]
        Pjν = [(-ν).exp()*sum([ν**(l*N+j)/factorial(l*N+j) for l in range(number_of_states)]
                            )for j in range(N)]
        for j in range(N):

            # μν_pow = (μ * ν).sqrt()  # Square root of the product of μ and ν
            FjZμXν[j] = (-(μ + ν) / 2).exp() / (Pjμ[j] * Pjν[j]).sqrt() * sum(
                [((μ * ν) ** Decimal((l * N + j) / 2) / factorial(l * N + j) ** Decimal('1.5'))
                    for l in range(number_of_states)])
            FjZμXν[j] = Decimal('1')-FjZμXν[j]
            FjXμXν[j] = (-(μ + ν) / Decimal('2')).exp() / (Pjμ[j] * Pjν[j]).sqrt() * sum(
                [((μ * ν) ** (Decimal(l * N + j)/Decimal('2')) / factorial(l * N + j))
                    for l in range(number_of_states)])
            FjXμXν[j] = Decimal('1')-FjXμXν[j]
        Y.append(FjXμXν[0])
plt.semilogy(range(2,32,2), Y)
plt.show()
