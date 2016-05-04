from numba import jitclass #import decorator
import numpy as np
from scipy.special import gammaln, digamma
from scipy.optimize import minimize


class simple_adaptive_rejection_sampler(object):
    """
    Parameters::
    @ lower_bound :  float 32
    @ upper_bound : float 32
    @ X : array of floats
    Description: 
    This is a template for an adapative rejection sampler implemented in Numba.
    Numba compiles python to machine code for significant speedups over Python.
    The only downside is that it's still early days for this module.
    In particular, current Numba versions do not accept functions as inputs
    This means that the function you want to sample from must be hard coded
    to the class.
    
    In particular, you must change these definitions for your purposes
    
    @ h : the log pdf of what you want to sample from, up to a normalizing constant 
    @ hprime: 
    """
    def __init__(self,lower_bound,upper_bound,X):
        """
        Immutable features 
        """
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.X = np.sort(X)
        self.b = X.size
        self.loga = np.sum(np.log(X))
        
    def h(self,alpha):
        return (alpha-1)*self.loga - self.b*gammaln(alpha)
        
    def hprime(self,alpha):
        return self.loga - self.b*digamma(alpha)
     
    @property
    def T(self):
        """
        These are the absicae that we are searching voer 
        """
        #this will define our envelope functions
        out = [self.X[0],self.X[self.b/2]]
        outer = self.X[self.b-1]
        while self.hprime(outer)>=0:
            outer = 2*outer
        out.append(outer)
        return out
    
    @property
    def hT(self):
        return [self.h(i) for i in self.T]
    
    @property
    def hprimeT(self):
        return [self.hprime(i) for i in self.T]
    
    @property
    def Z(self):
        """
        These define the critical points
        """
        for j in xrange(0,2):
            out.append(
            (self.hT[j+1] - self.hT[j] - \
                self.T[j+1]*self.hprimeT[j+1] + self.T[j]*self.hprimeT[j]) / 
            (self.hprimeT[j]-self.hprimeT[j+1]))
        return out
    
    def piecewise_upper(self,v):
        #determine which interval we're in, 0th, 1st, 2nd, etc.
        #the 0th would go from the lower bound to the first tangent intersection
        #the 1st would go from first intersection to the second etc.
        if v < self.Z[0]:
            return self.hT[0] + (v-self.T[0])*self.hprimeT[0]
        else:
            for i in self.Z):
        
        for i in range(len(self.Z)-1,-1,-1):
            if v >= self.Z[i]:
                break
            else: 
                continue   
        #we're in the ith, therefore
        
        
    def piecewise_lower(self,v):
        #determine which interval we're in, 0th, 1st, 2nd, etc.
        #the 1st would go from first intersection to the second etc.
        if v < self.T[0]:
            return -np.inf
        elif v > self.T[-1]:
            return -np.inf
        else:
            for i in xrange(0,len(self.T)-1):
                if v>=self.T[i+1]:
                    continue
                else:
                    break
            return ((self.T[i+1]-v)*self.hT[i]+
                    (v-self.T[i])*self.hT[i+1])/(self.T[i+1] - self.T[i])
        #we're in the ith, therefore
        

np.random.seed(1)
X = np.random.gamma(10,size=100)
lower_bound = 0.
upper_bound = np.inf
instance= simple_adaptive_rejection_sampler(lower_bound,upper_bound,X)
self=instance

x  = np.arange(0,30,0.25)
y = instance.h(x)
import matplotlib.pyplot as plt
plt.close("all")
plt.plot(x,y,instance.T,instance.hT,'bs',\
    x,[instance.piecewise_upper(i) for i in x],"r",
    x,[instance.piecewise_lower(i) for i in x],"y")
plt.show()
#print instance.D

plt.plot(x,y,instance.T,instance.hT,'bs',\
    x,[instance.piecewise_upper(i) for i in x],"r")

#instance.sample()
#print instance.D
#instance.sample()
##
#
#    def h(self,alpha):
#        return (alpha-1)*self.loga - self.b*np.math.lgamma(alpha)
#    
#    @property
#    def D(self):
#        return [self.X[0],self.X[self.b/2],self.X[self.b-1]]
#      
#    @D.setter  
#    def changeD(self,value):
#        self.D = 
#
#        

#
#instance= adaptive_rejection_sampler_template(lower_bound,upper_bound,X)
#print instance.X, instance.loga, instance.bi
#print instance.h(1),instance.h(2),instance.h(10)