# Adapative Rejection Sampler V 0.01

This is a Python implementation of the algorithm described in P. Wild, W.R. Gilks, Algorithm AS 287: Adaptive Rejection Sampling from Log Concave Density functions. Assuming that you have a log concave density, this algorithm allows you to sample without ever having to calculate an integral.  For the Gibbs sampling implementation this was created for, this is particularly attractive.

Below are kernel density plots from 10000 samples of the standard normal, Gamma (9,0.5), and truncated normal distributions.  I've also added used scipy to plot the pdf of these distributions.  As one can see, it's hard to tell which is which, so the sampler does appear to work.  

![alt text](https://github.com/libgober/ARS/blob/master/KDEplots1.png "KDE plots")

I calculated p-values for the ks-test of each set of samples against the hypothesized distribution. The smallest p-value of these was around 0.86.


In the bottom right quadrant I've done a small simulation study.  I drew 750 samples from a gamma distribution with shape 10 and scale 1, and then treat that as data for trying to perform Bayesian inference.  To do that, I've drawn 1000 samples from the posterior distribution of the conjugate prior to a gamma with unknown shape, scale = 1, and hyperpriors chosen so that the simulated data totally dominate. Implicitly, the distribution involves a parameter that is the product of all our data (in this example should be around 10^750). Using adaptive rejection sampling one can hit overflow very fast if one is not careful!  Most of the work one does with adjaptive rejection sampling happens in the log space, however, so if one tinkers with the algorithm a tiny bit one can mitigate such problems.  Indeed, this sampler was coded with these kind of problems in mind and seems to address this problem succesfully. As can be seen, the mode of the posterior distribution is close to its true value, 10, and all the mass is between 9.8 and 10.3. The Bayesian inference seems to have worked on a distribution whose parameter is treated by my laptop as inf!

The impetus for coding this sampler was to do fast Gibbs sampling on this very distribution, which has no closed form integral (the conjugate prior of a gamma likelihood with unknown shape). The concerns of this particular project dominated development, but future projects may again require such a sampler and I may do more work on it. Because Gibbs sampling need a single sample at a time, there did not seem to be much advantage in developing the lower hull described in the paper.  This is a feature that may be added at some point and it would presumably be fairly straightforward.  Hopefully the numba package will continue to make improvements and one day we'll get a big free speedup.  This module is mostly compliant with numba, some part of it is not (I forget which). At some point it may be translated in whole or in part to Cython, again as the project warrants.

Because the Bayesian posterior we are working with is very strongly peaked, the code presented here contains some modest enhancements to improve numerical stability in the case of overflow.  These are pretty straightforward, the basic idea is that wherever the paper suggests something like $exp(h(x))$, instead you try to do $exp(h(x)-h^*)$ where $h^*$ is the supremum of the upper envelope. The wisedom of this is apparent since you don't want to do $np.exp(1000) - np.exp(999)$. Working through this is not too difficult. 

Documentation is still pretty weak, but I'm not sure whether this is of interest to anyone else so I'll keep it that way for now.

