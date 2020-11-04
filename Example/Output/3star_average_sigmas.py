import sys
import numpy as np
from decimal import Decimal
import math
import pylab

INPUTNAME = "new-mcmc_fit.dat"
OUTPUTNAME = 'vals.dat'

input = open(INPUTNAME)
outputString = ""

for line in input:
    args = line.split(  )
    outputString+="{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18}\n".format(args[0], args[1], args[2], args[3], args[4], args[5],
            np.absolute(float(args[0])-float(args[2])), np.absolute(float(args[1])-float(args[3])), np.absolute(float(args[0])-float(args[4])), np.absolute(float(args[1])-float(args[5])),
            args[6], args[7], args[8], args[9], args[10], args[11], args[12], args[13], args[14])
        
with open (OUTPUTNAME, 'a') as f: f.write (outputString)


x1 = []
y1 = []
x2 = []
y2 = []
x3 = []
y3 = []
x1x2 = []
y1y2 = []
x1x3 = []
y1y3 = []
sep12 = []
sep13 = []
f1 = []
f2 = []
z = []
z1 = []
z2 = []
z3 = []
chi2 = []

File = open('vals.dat')
outputString=""

i=0
cuts=0
for plotPair in File:
    xandy = plotPair.split()
    x1.append(float(xandy[0]))
    y1.append(float(xandy[1]))
    x2.append(float(xandy[2]))
    y2.append(float(xandy[3]))
    x3.append(float(xandy[4]))
    y3.append(float(xandy[5]))
    x1x2.append(float(xandy[6]))
    y1y2.append(float(xandy[7]))
    x1x3.append(float(xandy[8]))
    y1y3.append(float(xandy[9]))
    sep12.append(float(xandy[10]))
    sep13.append(float(xandy[11]))
    f1.append(float(xandy[12]))
    f2.append(float(xandy[13]))
    z.append(float(xandy[14]))
    z1.append(float(xandy[15]))
    z2.append(float(xandy[16]))
    z3.append(float(xandy[17]))
    chi2.append(float(xandy[18]))

x1_avg = np.average(x1)
x1_std = np.std(x1)
y1_avg = np.average(y1)
y1_std = np.std(y1)
x2_avg = np.average(x2)
x2_std = np.std(x2)
y2_avg = np.average(y2)
y2_std = np.std(y2)
x3_avg = np.average(x3)
x3_std = np.std(x3)
y3_avg = np.average(y3)
y3_std = np.std(y3)
x1x2_avg = np.average(x1x2)
x1x2_std = np.std(x1x2)
y1y2_avg = np.average(y1y2)
y1y2_std = np.std(y1y2)
sep12_avg = np.average(sep12)
sep12_std = np.std(sep12)
sep13_avg = np.average(sep13)
sep13_std = np.std(sep13)
f1_avg = np.average(f1)
f1_std = np.std(f1)
f2_avg = np.average(f2)
f2_std = np.std(f2)
z_avg = np.average(z)
z_std = np.std(z)
z1_avg = np.average(z1)
z1_std = np.std(z1)
z2_avg = np.average(z2)
z2_std = np.std(z2)
z3_avg = np.average(z3)
z3_std = np.std(z3)
chi2_avg = np.average(chi2)
chi2_std = np.std(chi2)

print("x1 ",x1_avg,x1_std)
print("y1 ",y1_avg,y1_std)
print("x2 ",x2_avg,x2_std)
print("y2 ",y2_avg,y2_std)
print("x3 ",x3_avg,x3_std)
print("y3 ",y3_avg,y3_std)
print("1-2sep ",sep12_avg,sep12_std)
print("1-3sep ",sep13_avg,sep13_std)
print("f1 ",f1_avg,f1_std)
print("f2 ",(1.0-f1_avg),f1_std)
print("f3 ",(1.0-f1_avg),f1_std)
print("FTOTAL ",z_avg,z_std)
print("F1 ",z1_avg,z1_std)
print("F2 ",z2_avg,z2_std)
print("F3 ",z3_avg,z3_std)
print("chi2 ",chi2_avg,chi2_std)
#print("x1-x2 ",9.94*x1x2_avg,9.94*x1x2_std) #9.94=Keck pixel scale
#print("y1-y2 ",9.94*y1y2_avg,9.94*y1y2_std)
