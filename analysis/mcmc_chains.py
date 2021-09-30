##################
# Calculate average and standard deviations
# for parameters from mcmc_fit.dat output.
# Plot and save each chain.
##################
import os
import sys
import math
import pdb
import numpy as np
import matplotlib.pyplot as plt


def twostar(dir_mcmc, dir_out):
	print("\nStart Analysis")
	print('----Plotting chains')
	data = np.genfromtxt(dir_mcmc)
	params = ['X1', 'Y1', 'X2', 'Y2', 'Separation', 'FluxRatio',
	'FTOTAL', 'F1', 'F2', 'CHISQ']
	for i in range(np.shape(data)[1]):
		fig = plt.figure(figsize=(8, 6))
		ax = fig.add_subplot(111)

		plt.plot(np.arange(0,len(data)), data[:,i])
		plt.xlabel('step', fontsize=20)
		plt.ylabel(params[i], fontsize=20)
		plt.tick_params(which='both', length=8, width=1, labelsize=20, direction='in', right=True, top=True)

		plt.tight_layout()
		plt.savefig(dir_out + params[i] + '_chain.pdf')
		plt.close(fig)


def threestar(dir_mcmc, dir_out):
	print("\nStart Analysis")
	print('----Plotting chains')
	data = np.genfromtxt(dir_mcmc)
	params = ['X1', 'Y1', 'X2', 'Y2', 'X3', 'Y3', 'Separation',
	'Flux1', 'Flux2', 'FTOTAL', 'F1', 'F2', 'F3', 'CHISQ']
	for i in range(np.shape(data)[1]):
		fig = plt.figure(figsize=(8, 6))
		ax = fig.add_subplot(111)

		plt.plot(np.arange(0,len(data)), data[:,i])
		plt.xlabel('step', fontsize=20)
		plt.ylabel(params[i], fontsize=20)
		plt.tick_params(which='both', length=8, width=1, labelsize=20, direction='in', right=True, top=True)

		plt.tight_layout()
		plt.savefig(dir_out + params[i] + '_chain.pdf')
		plt.close(fig)


def avg_err_2star(dir_mcmc, dir_out, N_burn):
	print('----Calculating averages and errors')
	data = np.genfromtxt(dir_mcmc)
	params = ['X1', 'Y1', 'X2', 'Y2', 'Separation[mas]', 'FluxRatio',
	'FTOTAL', 'F1', 'F2', 'CHISQ']
	avgs = []
	errs = []
	xsep_avg = []
	ysep_avg = []
	xsep_err = []
	ysep_err = []
	for i in range(np.shape(data)[1]):
		avgs.append(np.average(data[N_burn:len(data),i]))
		errs.append(np.std(data[N_burn:len(data),i]))
	avgs_errs = np.stack((params, avgs, errs), axis=1)
	#--------------------------------------------------------
	for i in range(len(data)):
		xsep_avg.append((data[i,2]-data[i,0]))
		ysep_avg.append((data[i,3]-data[i,1]))
		xsep_err.append((data[i,2]-data[i,0]))
		ysep_err.append((data[i,3]-data[i,1]))
	xsep_avg = np.average(np.asarray(xsep_avg[N_burn:len(xsep_avg)]))*9.942
	ysep_avg = np.average(np.asarray(ysep_avg[N_burn:len(xsep_err)]))*9.942
	xsep_err = np.std(np.asarray(xsep_err[N_burn:len(xsep_err)]))*9.942
	ysep_err = np.std(np.asarray(ysep_err[N_burn:len(ysep_err)]))*9.942
	xyseps = np.vstack((xsep_avg, ysep_avg))
	xyseperrs = np.vstack((xsep_err, ysep_err))
	params2 = ('X2-X1', 'Y2-Y1')
	params2 = np.asarray(params2).reshape(2,1)
	#params3 = params2.reshape(2,1)
	xyseps = list(xyseps)
	xyseperrs = list(xyseperrs)
	final_seps = np.stack((params2, xyseps, xyseperrs), axis=1)
	#np.savetxt('./output/seps.dat', finallist)
	print(avgs_errs.shape)
	print(final_seps.shape)
	#final_list = np.append(avgs_errs, final_seps, axis=0)


	#np.savetxt(dir_out + 'avg_err_' + str(N_burn) + 'burn.dat', avgs_sigs, fmt='%s')

def avg_err_3star(dir_mcmc, dir_out, N_burn):
	print('----Calculating averages and errors')
	data = np.genfromtxt(dir_mcmc)
	params = ['X1', 'Y1', 'X2', 'Y2', 'X3', 'Y3', 'Separation',
	'Flux1', 'Flux2', 'FTOTAL', 'F1', 'F2', 'F3', 'CHISQ']
	avgs = []
	sigs = []
	for i in range(np.shape(data)[1]):
		avgs.append(np.average(data[N_burn:len(data),i]))
		sigs.append(np.std(data[N_burn:len(data),i]))
	avgs_sigs = np.stack((params, avgs, sigs), axis=1)
	np.savetxt(dir_out + 'avg_err_' + str(N_burn) + 'burn.dat', avgs_sigs, fmt='%s')






