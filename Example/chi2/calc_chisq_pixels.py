import os
import sys
import numpy as np
from decimal import Decimal
import math
import matplotlib.pyplot as plt

chi2s = []
xychi2 = []
chi2_rad = []
vals = []

box_size = 20 #Sqrt of fitting box area (rounded up to nearest digit).
x_center = [] #Central X position (i.e. between source/lens x-positions)
y_center = [] #Central Y position (i.e. between source/lens y-positions)

data = np.loadtxt('chisq_pixel.dat')

for i in range(1,250):
    for j in data:
        if abs(j[0]-1112) <= i and abs(j[1]-1073) <= i:
            chi2s.append(j)

chi2s = np.vstack(chi2s)
xychi2 = np.unique(chi2s, axis=0)

for m in range(2,box_size):
	if m % 2 != 0:
		vals.append([xychi2[0:(m**2),2].sum(),(m**2),xychi2[0:(m**2),2].sum()/(m**2)])

np.savetxt('chisq_vals.dat', vals, delimiter=' ', fmt='%s', header='chisq   total_pix   chi2/pix')

print(np.sum(xychi2[:,2]))


vals = np.asarray(vals)

fig = plt.subplots(figsize=(12,10))
ax1 = plt.subplot(111)

ax1.plot(np.sqrt(vals[:,1]), vals[:,0], 'k')

ax1.set_xlabel('Radius [pix]', fontsize=24)
ax1.set_ylabel(r'Total $\chi^2$', fontsize=24)
plt.tick_params(which='both', length=8, width=1, labelsize=22, direction='in', right=True, top=True)
plt.savefig('chisq_total.pdf')
#plt.show()

fig = plt.subplots(figsize=(12,10))
ax2 = plt.subplot(111)

ax2.plot(np.sqrt(vals[:,1]), vals[:,2], 'k')
plt.ylim(0.0,0.5)

ax2.set_xlabel('Radius [pix]', fontsize=24)
ax2.set_ylabel(r'$\chi^2$/pix', fontsize=24)
plt.tick_params(which='both', length=8, width=1, labelsize=22, direction='in', right=True, top=True)
plt.savefig('chisq_per_pix.pdf')
#plt.show()



