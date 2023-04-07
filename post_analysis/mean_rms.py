import sys
import os
import numpy as np
from decimal import Decimal
import math
import pylab
import pandas as pd
import statistics as st
from uncertainties import ufloat
from uncertainties.umath import *
import fnmatch
import pdb
import matplotlib.pyplot as plt
import matplotlib
from tabulate import tabulate

#Reading in the raw mcmc files

path = os.getcwd()
file_list = [k for k in os.listdir(path) if k.endswith('.mcmc')]
#pdb.set_trace()

if True:
    names = "X1_CENTER Y1_CENTER X2_CENTER Y2_CENTER X3_CENTER Y3_CENTER S1-2_SEP S1-3_SEP S1_F_CONTRIB S3_F_CONTRIB FTOTAL F1 F2 F3 CHISQ".split()
    use_columns = range(15)
    dtype = dict()
    [dtype.update({names[a]: np.float64}) for a in use_columns]
    column_name = [names[a] for a in use_columns]

    for i in range(len(file_list)):
        if i==0:
            output_mcmc = pd.read_csv(file_list[i], sep='\s+', header=None, names=column_name, 
                usecols=use_columns, dtype=dtype, skiprows=1)
            output_mcmc['ID'] = i
        else:
            data = pd.read_csv(file_list[i], sep='\s+', header=None, names=column_name, 
                usecols=use_columns, dtype=dtype, skiprows=1)
            data['ID'] = i
            output_mcmc = pd.concat([output_mcmc, data], ignore_index=True)

    output_mcmc.to_hdf('data.h5', key='df', mode='w')

else:
    output_mcmc = pd.read_hdf('data.h5')
    print(output_mcmc.head(4))

#Find the mean and rms of the deviation from the mean. 
#This is using the 'raw' mcmc files. Need to have a 'burn' value to remove from the beginning of the mcmc when it hasn't converged yet. 
#Dave recommends removing lines up until the chi^2 reaches within 1-2 of the chi^2 minimum. For me thats about 2.2%

def mean_rms(i, burn):
    mask_i = output_mcmc['ID'] == i
    #burn_mask_i = output_mcmc[250:len(mask_i)]
    
    DX21 = output_mcmc.loc[mask_i, 'X2_CENTER'] - output_mcmc.loc[mask_i, 'X1_CENTER']
    mean_dx21 = np.mean(DX21[burn:len(DX21)])
    dev_dx = DX21[burn:len(DX21)] - mean_dx21
    rms_dx = np.sqrt(np.mean(dev_dx**2))
    #print("mcmc", i, ": mean DX21 = ", np.mean(DX21), "and RMS = ", rms_dx)
    
    DY21 = output_mcmc.loc[mask_i, 'Y2_CENTER'] - output_mcmc.loc[mask_i, 'Y1_CENTER']
    mean_dy21 = np.mean(DY21[burn:len(DY21)])
    dev_dy = DY21[burn:len(DY21)] - mean_dy21
    rms_dy = np.sqrt(np.mean(dev_dy**2))
    #print("mcmc", i, ": mean DY21 = ", np.mean(DY21), "and RMS = ", rms_dy)
    
    Dr21 = np.sqrt((DX21**2)+(DY21**2))
    mean_dr21 = np.mean(Dr21[burn:len(Dr21)])
    dev_dr = Dr21[burn:len(Dr21)] - mean_dr21
    rms_dr = np.sqrt(np.mean(dev_dr**2))
    #print("mcmc", i, ": mean Dr21 = ", np.mean(Dr21), "and RMS = ", rms_dr)
    #print(dev)
    
    S2_F_CONTRIB = 1 - output_mcmc.loc[mask_i, 'S1_F_CONTRIB'] - output_mcmc.loc[mask_i, 'S3_F_CONTRIB']
    S21_FLUX_R = S2_F_CONTRIB / output_mcmc.loc[mask_i, 'S1_F_CONTRIB']
    mean_f21 = np.mean(S21_FLUX_R[burn:len(S21_FLUX_R)])
    dev_f12 = S21_FLUX_R[burn:len(S21_FLUX_R)] - mean_f21
    rms_f12 = np.sqrt(np.mean(dev_f12**2))
    #print("mcmc", i, ": mean f21 = ", np.mean(S21_FLUX_R), "and RMS = ", rms_f12)

    #plotting histogram
    #This applies to any of the parameters, just change "Dr21" part to whatever
    plt.figure(figsize=(10,6))
    hfont = {'fontname':'Helvetica'}
    matplotlib.rc('ytick', labelsize=15) 
    matplotlib.rc('xtick', labelsize=15) 
    hist = plt.hist(Dr21, bins='auto', color='#5070ff', rwidth=1)
    #plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Dr21 (pixels)', fontsize = 20, **hfont) 
    plt.ylabel('Frequency', fontsize = 20, **hfont)
    plt.savefig('Dr12.png', dpi=600)

    plt.clf()
    #plotting relation between (any) two parameters
    hfont = {'fontname':'Helvetica'}
    plt.figure(figsize=(10,8))
    plt.plot(DX21, DY21, '.', color = '#0a465d', markersize = 0.3)
    plt.xlabel('X2-X1', fontsize = 20, **hfont)
    plt.ylabel('Y2-Y1', fontsize = 20, **hfont)
    #plt.xlim(0,5)
    #plt.ylim(0,)
    #plt.axis('equal')
    plt.savefig('DX21_v_DY21.png', dpi=600)

    return round(mean_dx21,6), round(rms_dx,6), round(mean_dy21,6), round(rms_dy,6), round(mean_dr21,6), round(rms_dr,6), round(mean_f21,6), round(rms_f12,6)

 
# ~~~~~~~~~~~       
#Saving the function output (could probably put this in the function...)
dataframe=[]
for i in range(0,1): # range is the number of jackknife frames you have
    dataframe.append(mean_rms(i, 2100))
    x = dataframe
#print(x)
c = ['mean_dx','rms_dx','mean_dy21','rms_dy','mean_dr21','rms_dr21','mean_f_ratio','rms_fratio']
df=pd.DataFrame(x, columns=c)
dataframe.append(df)

df.to_csv('1195_mean_rms_3star_burn2100.dat',sep = '\t', header=True, mode='w')


# ~~~~~~~~~~

#Best Fit Values

path = os.getcwd()
file_list_bf = [k for k in os.listdir(path) if k.endswith('.bf')]
print(file_list_bf)
print(len(file_list_bf))

#Same as before (thank you Clement!)
if True:
    names = "X1_CENTER Y1_CENTER X2_CENTER Y2_CENTER X3_CENTER Y3_CENTER S1-2_SEP S1-3_SEP S1_F_CONTRIB S3_F_CONTRIB FTOTAL F1 F2 F3 CHISQ".split()
    use_columns = range(15)
    dtype = dict()
    [dtype.update({names[a]: np.float64}) for a in use_columns]
    column_name = [names[a] for a in use_columns]

    for i in range(len(file_list_bf)):
        if i==0:
            output_mcmc_bf = pd.read_csv(file_list_bf[i], sep='\s+', header=None, names=column_name, 
                usecols=use_columns, dtype=dtype, skiprows=1)
            output_mcmc_bf['ID'] = i
        else:
            data_bf = pd.read_csv(file_list_bf[i], sep='\s+', header=None, names=column_name, 
                usecols=use_columns, dtype=dtype, skiprows=1)
            data_bf['ID'] = i
            output_mcmc_bf = pd.concat([output_mcmc_bf, data_bf], ignore_index=True)

    #output_mcmc_bf['S2_F_CONTRIB'] = 1 - output_mcmc_bf['S1_F_CONTRIB'] - output_mcmc_bf['S3_F_CONTRIB']

    output_mcmc_bf.to_csv('bf_data.dat', sep = '\t', header=True, mode='w')

else:
    output_mcmc_bf = pd.read_csv('bf_data.dat')
    print(output_mcmc_bf.head(4))


#This takes the mcmc best-fit values, finds their mean, and then finds the jackknife error. I think for the errors the above function 
#(mean_rms()) is better?

N=1
#First get the average
seps_bf = np.array([np.mean(output_mcmc_bf.loc[output_mcmc_bf['ID'] == i, 'S1-2_SEP']) 
    for i in range(len(file_list_bf))])
avgs_sep_bf = np.mean(seps_bf)
print(avgs_sep_bf)
sep_jkkf = np.sqrt(((N-1)/N)*sum((seps_bf-np.mean(seps_bf))**2))
print("average S21 separation=",avgs_sep_bf, "+/-", sep_jkkf, "mas")                 
#np.sqrt(((N-1)/N)*sum((std_DX-np.mean(std_DX))**2))


DX1_bf = np.array([np.mean(output_mcmc_bf.loc[output_mcmc_bf['ID'] == i, 'X1_CENTER']) 
    for i in range(len(file_list_bf))])
DX2_bf = np.array([np.mean(output_mcmc_bf.loc[output_mcmc_bf['ID'] == i, 'X2_CENTER']) 
    for i in range(len(file_list_bf))])
DX_bf=DX1_bf-DX2_bf
avgs_DX = np.mean(DX_bf)
DX_jkkf = np.sqrt(((N-1)/N)*sum((DX_bf-np.mean(DX_bf))**2))
print("average S21 DX=", avgs_DX*9.942, "+/-", DX_jkkf*9.942, "mas")

DY1_bf = np.array([np.mean(output_mcmc_bf.loc[output_mcmc_bf['ID'] == i, 'Y1_CENTER']) 
    for i in range(len(file_list_bf))])
DY2_bf = np.array([np.mean(output_mcmc_bf.loc[output_mcmc_bf['ID'] == i, 'Y2_CENTER']) 
    for i in range(len(file_list_bf))])
DY_bf=DY1_bf-DY2_bf
avgs_DY = np.mean(DY_bf)
DY_jkkf = np.sqrt(((N-1)/N)*sum((DY_bf-np.mean(DY_bf))**2))
print("average S21 DY=", avgs_DY*9.942,"+/-", DY_jkkf*9.942, "mas")

F1_bf = np.array([np.mean(output_mcmc_bf.loc[output_mcmc_bf['ID'] == i, 'F1']) 
    for i in range(len(file_list_bf))])
F2_bf = np.array([np.mean(output_mcmc_bf.loc[output_mcmc_bf['ID'] == i, 'F2']) 
    for i in range(len(file_list_bf))])
FLX_bf=F2_bf/F1_bf
avgs_FLX = np.mean(FLX_bf)
FLX_jkkf = np.sqrt(((N-1)/N)*sum((FLX_bf-np.mean(FLX_bf))**2))
print("average S2/S1 Flux Ratio=", avgs_FLX, "+/-", FLX_jkkf)



