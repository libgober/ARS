# -*- coding: utf-8 -*-
"""
Created on Wed May  4 18:47:16 2016

@author: brianlibgober
"""
#%%
from ars_libgober import *
#%%

np.random.seed(1)
X = np.sort(np.random.gamma(10,size=300))
loga = np.log(X).sum(); b = X.size
from scipy.special import digamma
def h(alpha):
    return (alpha-1)*loga - b*np.math.lgamma(alpha)
h = np.vectorize(h)
def hprime(alpha):
    return loga - b*digamma(alpha)
hprime = np.vectorize(hprime)
initial_knots = np.array([X[0],X[50],X[-1]])
xlb = 0
xub = np.inf
np.random.seed(1)
#%% setup
example = ArsSampler(initial_knots,np.vectorize(h),np.vectorize(hprime),0,xub)
example.sample(1)