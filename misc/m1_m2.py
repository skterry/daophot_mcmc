#============================================================
# Given a flux ratio and total magnitude for two stars,
# Calculate the individual magnitudes for each star.
# Run script from the terminal
#============================================================
import sys
import numpy as np

mtot   = float(input("Enter Total Magnitude: "))
fratio = float(input("Enter Flux Ratio (F1/F2): "))

F2 = (10**(-0.4*mtot))/(fratio+1)
F1 = fratio*F2

m1 = -2.5*np.log10(F1)
m2 = -2.5*np.log10(F2)

print("-----------------------------------")
print("m1 =", m1)
print("m2 =", m2)