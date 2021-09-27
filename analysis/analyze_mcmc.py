"""
################################
Analyze output from DAOPHOT-MCMC
################################
Parameters:
dir_mcmc (string) - location of 'mcmc_fit.dat' file
dir_out (string) - location of output values and plots 
"""

from analysis.mcmc_chains import twostar, threestar

dir_mcmc = '/Users/seanterry/keck-code/ob110950_2019/combo_a_2021_09_24/mcmc_fit.dat'
dir_out = './output/'

twostar(dir_mcmc, dir_out)