# phylo-particleMCMC
Matlab code for phylodynamic inference using particle filtering and particle MCMC. Implemented in Matlab. See Rasmussen et al. (2011) for a full description of the algorithms: https://doi.org/10.1371/journal.pcbi.1002136

To run the particle MCMC algorithm, run the main_Inference() function in Matlab. main_Inference calls run_MCMC to run the MCMC. Each MCMC step, get_Likelihoods is called to compute the marginal likelihood which in turn calls run_SMC to run the particle filtering algorithm.

The data files SIR_endemic_sim3_mockData.mat and SIR_endemic_sim3_coalTimes are provided as examples. These data files contain the mock time series and genealogy shown in Figures 1 and 5 in the paper.

Note to Octave users: While the code was never extensively tested in Octave, it should run with one exception. The Matlab function mvnrnd.m for multivariate normal random number generation is not supported in Octave but is available at: 

http://homepages.inf.ed.ac.uk/imurray2/code/matlab_octave_missing/mvnrnd.m
