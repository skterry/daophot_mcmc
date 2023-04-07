# DAOPHOT-MCMC

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

DAOPHOT-MCMC is a modified version of the subroutine NSTAR.F which implements a Markov chain Monte Carlo (MCMC) routine for fitting highly blended stellar positions and fluxes. This code is installed and implemented within the overall [DAOPHOT-II](https://www.star.bris.ac.uk/~mbt/daophot/) structure. Instructions and descriptions of the original software suite can be found in the DAOPHOT-II [User's Manual (direct download)](https://www.star.bris.ac.uk/~mbt/daophot/mud9.ps).

Some of the DAOPHOT-MCMC functionality includes:

* Simultaneous 1-star, 2-star, and 3-star PSF fitting.
* Metropolis-Hastings routine to calculate posterior probability distributions for positions, separations, flux ratios, and total flux.
* Best-fit ![equation](https://latex.codecogs.com/gif.latex?%5Cchi%5E2)/pixel map to measure PSF stability.
* (Optional) constraints on the total flux and separation(s) (source code edits required).
* (Optional) shell scripts (.sh) to automate MCMC runs on large datasets.

Detailed information can be found in the DAOPHOT-MCMC User's Manual [here](https://github.com/skterry/DAOPHOT-MCMC/blob/master/User_Manual.pdf).

If you find this code useful in your research, please cite [Terry et al. 2021](https://iopscience.iop.org/article/10.3847/1538-3881/abcc60)
