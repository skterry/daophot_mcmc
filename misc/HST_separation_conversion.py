#============================================================
# Given a flux ratio and total magnitude for two stars,
# Calculate the individual magnitudes for each star.
# Run script from the terminal
#============================================================
import sys
import numpy as np

theta_deg   = float(input("Enter HST rotation (hdr keyword is 'PA_APER'): "))
x1 = float(input("Enter X1 position: "))
y1 = float(input("Enter Y1 position: "))
x2 = float(input("Enter X2 position: "))
y2 = float(input("Enter Y2 position: "))

theta = np.deg2rad(theta_deg)

x1prime = x1*np.cos(theta) - y1*np.sin(theta)
y1prime = x1*np.sin(theta) + y1*np.cos(theta)

x2prime = x2*np.cos(theta) - y2*np.sin(theta)
y2prime = x2*np.sin(theta) + y2*np.cos(theta)

print("-----------------------------------")
print("New X1, Y1 =", x1prime, y1prime)
print("New X2, Y2 =", x2prime, y2prime)

print("-----------------------------------")
print("East separation (x1-x2) =", x1prime-x2prime)
print("North separation (y1-y2) =", y1prime-y2prime)

print("-----------------------------------")
print("East separation (mas) =", (x1prime-x2prime)*40)
print("North separation (mas) =", (y1prime-y2prime)*40)