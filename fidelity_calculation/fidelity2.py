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
    # with localcontext() as ctx:  # Use local context to avoid affecting global precision
    # ctx.prec = 100  # Set precision for this block
    μ = Decimal('0.5')
    # ν = Decimal('0.005')
    F = sum([sum([((-μ).exp()*μ**(l*N+j)/factorial(l*N+j))**Decimal('2') for l in range(number_of_states)]).sqrt() for j in range(N)])
    Y.append((Decimal(1)-F))
plt.figure(figsize=(6,4))
plt.semilogy(range(2,32,2), Y)
plt.xlabel('Number of Randomized Discrete Phases: N')
plt.ylabel('Difference between Eq.2 and Eq.3: (1-F^2)^0.5')
# plt.savefig('fig_fidelity2.eps',dpi=600)
plt.show()
