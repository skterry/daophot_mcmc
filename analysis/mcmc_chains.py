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
	ncols = np.shape(data)[1]
	params = ['X1', 'Y1', 'X2', 'Y2', 'Separation', 'FluxRatio',
	'FTOTAL', 'F1', 'F2', 'CHISQ']
	for i  in range(ncols):
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
	ncols = np.shape(data)[1]
	params = ['X1', 'Y1', 'X2', 'Y2', 'X3', 'Y3', 'Separation',
	'Flux1', 'Flux2', 'FTOTAL', 'F1', 'F2', 'F3', 'CHISQ']
	for i  in range(ncols):
		fig = plt.figure(figsize=(8, 6))
		ax = fig.add_subplot(111)

		plt.plot(np.arange(0,len(data)), data[:,i])
		plt.xlabel('step', fontsize=20)
		plt.ylabel(params[i], fontsize=20)
		plt.tick_params(which='both', length=8, width=1, labelsize=20, direction='in', right=True, top=True)

		plt.tight_layout()
		plt.savefig(dir_out + params[i] + '_chain.pdf')
		plt.close(fig)


def avg_sig_2star(dir_mcmc, dir_out):
	print('----Calculating averages and errors')
	data = np.genfromtxt(dir_mcmc)
	ncols = np.shape(data)[1]
	params = ['X1', 'Y1', 'X2', 'Y2', 'Separation', 'FluxRatio',
	'FTOTAL', 'F1', 'F2', 'CHISQ']
	avgs = []
	sigs = []

	for i  in range(ncols):
		avgs.append(np.average(data[:,i]))
		sigs.append(np.std(data[:,i]))
	avgs_sigs = np.stack((params, avgs, sigs), axis=1)
	np.savetxt(dir_out + 'avg_err.dat', avgs_sigs, fmt='%s')












