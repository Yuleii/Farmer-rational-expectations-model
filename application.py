# Import Python Packages
import numpy as np
from function import solve_farmer_model
import seaborn as sns
import matplotlib.pyplot as plt


# Set parameter values.
beta = 0.95
alpha = 0.3
delta = 0.1
rho = 0.9

# Set the number of time periods to simulate.
T = 100

# first 10 periods given
inital_productivity_state = [0, -0.005, -0.009,
                             -0.013, -0.022, -0.021, -0.019, -0.011, -0.012, -0.003]

# Get simulation result.
s_t, c_t, k_t, y_t = solve_farmer_model(alpha, beta, delta, rho, T, inital_productivity_state)

# plot time series.
time_1 = np.array(list(range(T)))
time_2 = np.array(list(range(T + 1)))

sns.set_style("whitegrid")
fig, ax = plt.subplots(figsize=(T / 8, T / 16))

# Consumption levels, capital levels, and output levels against time.
ax.plot(time_1, s_t, label="productivity state", color="C7", alpha=0.7, linestyle='dashdot')
ax.plot(time_1, c_t, label="consumption", color="C1", alpha=0.7, linestyle='dotted')
ax.plot(time_2, k_t, label="capital", color="C2", alpha=0.7, linestyle='dashed')
ax.plot(time_1, y_t, label="output", color="C0", alpha=0.7)
ax.legend(loc="upper left", fontsize=T / 10)

plt.xlim((0, T))
plt.xlabel("Time", fontsize=T / 6)
plt.title("Productivity state, consumption, capital and output over time", fontsize=T / 5)
fig.savefig("figures/time_series_plot")
