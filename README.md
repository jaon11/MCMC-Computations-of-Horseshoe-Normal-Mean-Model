# MCMC Computations of Horseshoe Normal Mean Model
This is a project example of a horseshoe prior with normal mean model computation. Here is an introduction of the files in this repository. 


## MCMC Computations of Horseshoe Normal Mean Model.pdf:
Studying note of MCMC Computations of Horseshoe Normal Mean Model.


## Codes:
### hs_nm.R:  
main function including 5 approaches to compute the horseshoe normal mean model.
	 1. slice: stepping-out slice sampling dealing with tau and lambda(others: Gibbs sampling)
	 2. MwG: Metropolis within Gibbs sampling (MH: lambda & tau)
	 3. Gibbs: Gibbs sampling using Wand mixture
	 4. MSwG: Metropolis slice within Gibbs sampling (MH: lambda; MH: tau)
	 5. slice2: Damien et al.(1999) slice sampling dealing with tau and lambda (others: Gibbs sampling)

### hs_nm_analysis: 
simulation of hs_nm(), comparing with HS.normal.means() in horseshoe r package, including computation stability, efficiency, signal-noise classification and figure results.


## Folders: 
**example 1:** single-signal stage: sigma2, tau ACF, and shrinkage figures with signal=2,5,10 in 10 seeds: 2004-2013

**example 2:** multi-signal stages: sigma2,tau ACF, and shrinkage figures with signal=2,5,10 in 10 seeds: 2004-2013


## txt:
**compu S&E.txt:** tables of computation stability, efficiency and signal-noise classification with 200 seeds: 2004-2203.
