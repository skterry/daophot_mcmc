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
from numpy.random import uniform
from numpy import convolve, cumsum, histogram, linspace
import matplotlib.pyplot as plt


def pdf(dir_mcmc, dir_out, N_jackknife, prefix, N_burn):
	chains_files = []
	for ee in range(N_jackknife):
		lis = np.genfromtxt(dir_mcmc + prefix + str(ee+1) + '_ks.mcmc')
		chains_files.append(lis)

	all_pdfs = []
	for i in enumerate(chains_files):
		pdf_sep= histogram(i[1][:,4], bins= 50)[0]
		all_pdfs.append(pdf_sep)

	fig = plt.figure(figsize=(8,7))
	ax1 = fig.add_subplot(111)
	for i in enumerate(all_pdfs):
		print(len(i[1]))
		#ax1.hist(i[1][:,0], bins=10)
		#plt.show()