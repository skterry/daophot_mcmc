"""
################################
Analyze output from DAOPHOT-MCMC (+Jackknife)
################################
Parameters:
dir_mcmc (string) - location of '*.mcmc' file(s)
dir_chisq_pix (string) - location of '*.chi2pix' file(s)
dir_out (string) - location of output values and plots

box_rad (int) - sqrt(AREA) of fitting box (box limits are printed to terminal at beginning of MCMC run)
i.e. if fit box X=[1108,1129], Y=[1055,1089], then box_rad = 27 (round up to nearest whole number)

N_burn (int) - Burn in amount
"""

from post_analysis.mcmc_chains import twostarchains, threestarchains, avg_err_2star, avg_err_3star
from post_analysis.convolve_pdfs import pdf
#from analysis.corner_plot import corner
#from analysis.jackknife_calc import calc

#Directories
dir_mcmc = '/Users/seanterry/keck-code/ob110950/ob110950_2019/a_2021_09_30_JackknifeFrames/'
dir_chisq_pix = '/Users/seanterry/keck-code/ob110950/ob110950_2019/a_2021_09_30_JackknifeFrames/'
dir_out = './output/'
prefix = 'magtest1'

#Parameters
N_jackknife = 21
box_rad = 20
N_burn = 250


twostarchains(dir_mcmc, dir_out, N_jackknife, prefix, N_burn)
#avg_err_2star(dir_mcmc, dir_out, N_jackknife, prefix, N_burn)
#chi2pix(dir_chisq_pix, dir_out, box_rad)

#pdf(dir_mcmc, dir_out, N_jackknife, prefix, N_burn)



#----------3 star MCMC-----------
#threestarchains(dir_mcmc, dir_out) #uncomment this for 3star
#avg_err_3star(dir_mcmc, dir_out) #uncomment this for 3star




