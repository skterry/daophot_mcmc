# DAOPHOT-MCMC

DAOPHOT-MCMC is a modified version of the subroutine NSTAR.F which implements a Markov chain Monte Carlo (MCMC) routine for fitting crowded stellar positions and fluxes. This code is installed and implemented within the overall [DAOPHOT-II](https://www.star.bris.ac.uk/~mbt/daophot/) structure (**note**: you must have a pre-installed version of DAOPHOT-II). Instructions and descriptions of the original software suite can be found in the DAOPHOT-II [User's Manual (direct download)](https://www.star.bris.ac.uk/~mbt/daophot/mud9.ps).

Some of the DAOPHOT-MCMC functionality includes:

* Simultaneous 1-star, 2-star, and 3-star PSF fitting.
* Metropolis-Hastings routine to deliver posterior probability distributions for positions, separations, flux ratios, and total flux.
* Best-fit $`\chi^2`$/pixel map.
* (Optional) constraints on the total flux and separation(s) for 2-star and 3-star fitting.
* (Optional) shell scripts (.sh) to automate MCMC runs on large datasets.

Detailed information can be found in the DAOPHOT-MCMC User's Manual [here](https://github.com/skterry/DAOPHOT-MCMC/blob/master/User_Manual.pdf).

A Python wrapper (with added functionality) has also been developed, you can find it [here](https://github.com/clementranc/nstarwrap).

If you find this code useful in your research, please cite [Terry et al. 2021](https://iopscience.iop.org/article/10.3847/1538-3881/abcc60)
