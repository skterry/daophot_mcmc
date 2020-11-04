import sys
import numpy as np
from decimal import Decimal
import math
import pylab

INPUTNAME = "mcmc_fit.dat"
OUTPUTNAME = 'temp-vals.dat'

input = open(INPUTNAME)
outputString = ""

for line in input:
    args = line.split(  )
    outputString+="{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n".format(args[0], args[1], args[2], args[3],
            np.absolute(float(args[0])-float(args[2])), np.absolute(float(args[1])-float(args[3])),
            args[4], args[5], args[6], args[7])
        
with open (OUTPUTNAME, 'a') as f: f.write (outputString)


x1 = []
y1 = []
x2 = []
y2 = []
x1x2 = []
y1y2 = []
sep = []
f1 = []
z = []
chi2 = []

File = open('temp-vals.dat')
outputString=""

i=0
cuts=0
for plotPair in File:
    xandy = plotPair.split()
    x1.append(float(xandy[0]))
    y1.append(float(xandy[1]))
    x2.append(float(xandy[2]))
    y2.append(float(xandy[3]))
    x1x2.append(float(xandy[4]))
    y1y2.append(float(xandy[5]))
    sep.append(float(xandy[6]))
    f1.append(float(xandy[7]))
    z.append(float(xandy[8]))
    chi2.append(float(xandy[9]))

x1_avg = np.average(x1)
x1_std = np.std(x1)
y1_avg = np.average(y1)
y1_std = np.std(y1)
x2_avg = np.average(x2)
x2_std = np.std(x2)
y2_avg = np.average(y2)
y2_std = np.std(y2)
x1x2_avg = np.average(x1x2)
x1x2_std = np.std(x1x2)
y1y2_avg = np.average(y1y2)
y1y2_std = np.std(y1y2)
sep_avg = np.average(sep)
sep_std = np.std(sep)
f1_avg = np.average(f1)
f1_std = np.std(f1)
z_avg = np.average(z)
z_std = np.std(z)
chi2_avg = np.average(chi2)
chi2_std = np.std(chi2)

print("x1 ",x1_avg,x1_std)
print("y1 ",y1_avg,y1_std)
print("x2 ",x2_avg,x2_std)
print("y2 ",y2_avg,y2_std)
print("sep ",sep_avg,sep_std)
print("f1 ",f1_avg,f1_std)
print("f2 ",(1.0-f1_avg),f1_std)
print("z ",z_avg,z_std)
print("chi2 ",chi2_avg,chi2_std)
print("x1-x2 ",9.94*x1x2_avg,9.94*x1x2_std) #9.94=Keck pixel scale
print("y1-y2 ",9.94*y1y2_avg,9.94*y1y2_std)
