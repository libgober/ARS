# Adapative Rejection Sampler V 0.01

This is a Python implementation of the algorithm described in P. Wild, W.R. Gilks, Algorithm AS 287: Adaptive Rejection Sampling from Log Concave Density functions. Assuming that you have a log concave density, this algorithm allows you to sample without ever having to calculate an integral.  For the Gibbs sampling implementation this was created for, this is particularly attractive.

The first three figures are kernel density plots from 10000 samples of the standard normal, Gamma (9,0.5), and truncated normal distributions.  I've also added used scipy to plot the pdf of these distributions.  As one can see, it's hard to tell which is which, so the sampler does appear to work.  

![alt text](https://github.com/libgober/ARS/blob/master/KDEplots.png "KDE plots")

I calculated p-values for the ks-test of each set of samples against the hypothesized distribution. The smallest p-value of these was around 0.86.

In the bottom right quadrant I've done a small simulation study.  I drew 100 samples from a gamma distribution with shape 10 and scale 1, and then treat that as data for trying to perform Bayesian inference.  To do that, I've drawn 1000 samples from the posterior distribution of the conjugate prior to a gamma with unknown shape, scale = 1, and hyperpriors chosen so that the simulated data totally dominate. Implicitly, the distribution involves a parameter that is the product of all our data (in this example should be around 10^100), so we can hit overflow very fast.  Since we are careful about our work in the log space, we don't have any problems. As can be seen, the mode of the posterior distribution is close to 10 and all the mass is between 9 and 11.  The Bayesian inference seems to have worked.

The impetus for coding this sampler was to do fast Gibbs sampling on a distribution that has no closed form integral (the conjugate prior of a gamma likelihood with unknown shape). The concerns of that project dominated development. Because we only need a single sample at a time, there did not seem to be much advantage in developing the lower hull described in the paper.  This is a feature that may be added at some point and it would presumably be fairly straightforward.  Hopefully the numba package will continue to make improvements and one day we'll get a big free speedup.  This module is mostly compliant with numba, some part of it is not (I forget which). At some point it may be translated in whole or in part to Cython.

Because the Bayesian posterior we are working with is very strongly peaked, the code presented here contains some modest enhancements to improve numerical stability in the case of overflow.  These are pretty straightforward, the basic idea is that wherever the paper suggests something like $exp(h(x))$, instead you try to do $exp(h(x)-h^*)$ where $h^*$ is the supremum of the upper envelope. The wisedom of this is apparent since you don't want to do $np.exp(1000) - np.exp(999)$. Working through this is not too difficult. 

Documentation is still pretty weak, but I'm not sure whether this is of interest to anyone else so I'll keep it that way for now.

