#============================================================
# Given a flux ratio and total magnitude for two stars,
# Calculate the individual magnitudes for each star.
# Run script from the terminal
#============================================================
import sys
import numpy as np

theta_deg = float(input("Enter Keck rotation (hdr keyword is 'ROTPOSN'): "))
dx12 = float(input("Enter X1-X2 Separation (in pix): "))
dy12 = float(input("Enter Y1-Y2 Separation (in pix): "))

theta = np.deg2rad(-1*theta_deg)

dx12prime = (dx12)*np.cos(theta) - (dy12)*np.sin(theta)
dy12prime = (dx12)*np.sin(theta) + (dy12)*np.cos(theta)

#dx2prime = x2*np.cos(theta) - y2*np.sin(theta)
#dy2prime = x2*np.sin(theta) + y2*np.cos(theta)

#print("-----------------------------------")
#print("New X1, Y1 =", x1prime, y1prime)
#print("New X2, Y2 =", x2prime, y2prime)

print("-----------------------------------")
print("East separation (x1-x2) =", dx12prime)
print("North separation (y1-y2) =", dy12prime)

print("-----------------------------------")
print("East separation (mas) =", (dx12prime)*9.942)
print("North separation (mas) =", (dy12prime)*9.942)