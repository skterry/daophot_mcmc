#============================================================
# Given a rotation angle (theta) and separation between two objects (in pixels),
# Calculate the separation of the two objects in East and North directions.
# Run script from the terminal
#============================================================
import sys
import numpy as np

theta_deg = float(input("Enter HST rotation (hdr keyword is 'PA_APER'): "))

dx12 = float(input("Enter X2-X1 Separation (in pix): "))
edx12 = float(input("Enter error (sigma) in X2-X1 Separation (in pix): "))

dy12 = float(input("Enter Y2-Y1 Separation (in pix): "))
edy12 = float(input("Enter error (sigma) Y2-Y1 Separation (in pix): "))

theta = np.deg2rad(theta_deg)

dx12prime = (dx12)*np.cos(theta) - (dy12)*np.sin(theta)
edx12prime = np.hypot(edx12*np.cos(theta),edy12*np.sin(theta))

dy12prime = (dx12)*np.sin(theta) + (dy12)*np.cos(theta)
edy12prime = np.hypot(edx12*np.sin(theta),edy12*np.cos(theta))

print("-----------------------------------")
print("East separation (mas) =", (dx12prime)*40, "+/-", (edx12prime)*40)
print("North separation (mas) =", (dy12prime)*40, "+/-", (edy12prime)*40)