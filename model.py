# MSEIRVSD Deterministic and Stochastic Model

import matplotlib.pyplot as plt
import numpy as np

# Defining the constants to be used
alpha = 1/5.8 # inverse of mean incubation period
beta = 0.3    # contact rate
gamma = 1/5   # inverse of mean recovery period
N = 70000000     # total population
CIR = 40      # inverse of cases to infections ratio
NDays = 600  # Number of days in the simulation

# Creating the NDays array
Days = np.linspace(1, NDays, NDays)

# Defining the compartment array
Compartment = np.zeros((6, NDays))

# Initialization

# Defining the seeding equation
for idx in range(NDays):
    if idx > 6:
        break
    Compartment[5][idx] += 100

# Defining the initial compartmental values
Compartment[0][0] = N - Compartment[5][0]
Compartment[1][0] = Compartment[5][0]

# The main loop
for idx in range(NDays-1):

    # Defining the dynamical equations
    Compartment[0][idx+1] = Compartment[0][idx] - beta * Compartment[0][idx] * (1/(N-1)) * Compartment[2][idx] - Compartment[5][idx]
    Compartment[1][idx+1] = Compartment[1][idx] - (Compartment[0][idx+1] - Compartment[0][idx]) - alpha * Compartment[1][idx] + Compartment[5][idx]
    Compartment[2][idx+1] = Compartment[2][idx] + alpha * Compartment[1][idx] - gamma * Compartment[2][idx]
    Compartment[3][idx+1] = Compartment[3][idx] + gamma * Compartment[2][idx]
    Compartment[4][idx] = (alpha/CIR) * Compartment[1][idx]
    Compartment[4][idx+1] = (alpha / CIR) * Compartment[1][idx+1]

# Answers to questions
print("beta_trial_and_error:", end=" ")
print(beta)
print("num_days_to_peak:", end=" ")
print(np.argmax(Compartment[4][:]))
print("num_days_greater_than_1000_cases:", end=" ")
print(sum(Compartment[4][:] > 1000))

# Plotting
plt.title("Karnataka COVID Epidemic")
plt.xlabel("Number of days elapsed since initial seeding")
plt.ylabel("Number of people in the compartments")
graphs = plt.plot(Days, np.rint(Compartment).transpose())
plt.legend(graphs, ("Susceptible", "Exposed", "Infectious", "Recovered", "Cases", "Seeds"))
plt.yscale('asinh')
# plt.yscale('linear')
plt.grid()
plt.show()