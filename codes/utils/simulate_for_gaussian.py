import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# Generate x values
x = np.linspace(-15, 0, 100)

# Calculate the CDF values for the standard normal distribution
cdf = norm.cdf(x)

# Plot the CDF using semilogy
plt.semilogy(x, cdf)

# Set the plot title and labels
plt.title('Standard Normal Distribution - Cumulative Distribution Function')
plt.xlabel('x')
plt.ylabel('CDF')

# Display the plot
plt.show()
