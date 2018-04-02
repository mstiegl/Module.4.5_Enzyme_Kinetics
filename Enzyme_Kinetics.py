import numpy as np
import matplotlib.pyplot as plt
# Simulation of enzyme kinetics

# Set up our time increment and our array of  time values
dt = .05   # s
t= np.arange(0,3,dt)

# Constants
k1 = 0.05       # 1/mM*s (rate of substrate removal per mM of enzyme)
k2 = 0.10       # 1/s (rate of substrate re-formation)
k3 = 0.02       # 1/s (rate of product formation)
k4 = 0.00       # 1/s (rate of reverse formation of intermediate complex)
initE = 1       # mM
initS = 0       # mM (we'll vary this along a range of values, below)
initES = 0      # mM
initP = 0       # mM


# Unlike most of the other examples, here we're going to run our simulation
# numerous times, using a different value of initS each time. This means we
# will have a nested-for loop: one for loop to iterate through the
# different starting values of S, and an inner for loop to go through the
# different time steps.

# The initSvals vector will contain the values we want to use for initS.
# The ratesOfP vector will contain our "answers" for each simulation. Every
# time we run the simulation with a value of initS, we'll record the result
# in ratesOfP at the appropriate spot (iterNum, which is an integer that
# increments one value for each simulation).
initSvals = np.arange(0,30,1.5,dtype=float)
ratesOfP = np.arange(len(initSvals),dtype=float)
iterNum = 0

for initS in initSvals:

    # Set up our stock variables and initial conditions.
    S = np.arange(len(t),dtype=float)
    S[0] = initS
    E =  np.arange(len(t),dtype=float)
    E[0] = initE
    ES =  np.arange(len(t),dtype=float)
    ES[0] = initES
    P = np.arange(len(t),dtype=float)
    P[0] = initP

    # Loop a standard number of times, starting with i=2 since we've already
    # set up the simulation's initial conditions at i=1.
    for i in range (1,len(t)):

        # Compute rates of change (in people per day) for each of the three
        #   categories. (We're just doing Euler's method here.)
        Sprime = k2 * ES[i-1] - k1 * E[i-1] * S[i-1]
        Eprime = k2 * ES[i-1] - k1 * E[i-1] * S[i-1] +k3 * ES[i-1] - k4 * E[i-1] * P[i-1]
        ESprime = k1 * E[i-1] * S[i-1] - k2 * ES[i-1] +k4 * E[i-1] * P[i-1] - k3 * ES[i-1]
        Pprime = k3 * ES[i-1] - k4 * E[i-1] * P[i-1]

        # Increase or decrease the population of each category based on
        #   the rates of change for this time step.
        S[i] = S[i-1] + (Sprime) * dt
        E[i] = E[i-1] + (Eprime) * dt
        ES[i] = ES[i-1] + (ESprime) *dt
        P[i] = P[i-1] + (Pprime) *dt
        print(P[i])


        # We're now done with one simulation. Add the final P' value to our
    # vector ratesOfP, which is collecting our results (one per simulation)
    # so we can plot them at the end.
        ratesOfP[iterNum]=Pprime*1000

    iterNum = iterNum + 1
print(ratesOfP)
plt.bar(initSvals,ratesOfP)
plt.show()

