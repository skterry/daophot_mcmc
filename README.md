# DAOPHOT-MCMC

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

DAOPHOT-MCMC is a modified version of the subroutine NSTAR.F which implements a Markov chain Monte Carlo (MCMC) routine for fitting highly blended stellar positions and fluxes. This code is installed and implemented within the overall [DAOPHOT-II](http://www.star.bris.ac.uk/~mbt/daophot/) structure. Instructions and descriptions of the original software suite can be found in the DAOPHOT-II [User's Manual](http://www.astro.wisc.edu/sirtf/daophot2.pdf).

Some of the DAOPHOT-MCMC functionality includes:

* Simultaneous 1-star, 2-star, and 3-star PSF fitting.
* Chi^2 minimization to determine best-fit position, separation, flux ratio, and total flux.
* Parameter constraints for best-fit MCMC.

The User's Manual for DAOPHOT-MCMC is a work in progress, and can be found [here](https://github.com/skterry/DAOPHOT-MCMC/blob/master/User_Manual.pdf).

An example MCMC simulation output can be seen below:

![gif 1](https://github.com/skterry/DAOPHOT-MCMC/blob/master/Example_output/gif-change.gif)

