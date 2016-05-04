from numba import jitclass #import decorator
from numba import float32,float64,int32,int64 #import types
import numpy as np
from scipy.special import digamma,gammaln



spec = [
('lower_bound', float32),
('upper_bound',float32),
('X', float64[:]),
('D',float64[:]),
('b',int32),
('loga',float64)
]

@jitclass(spec)
class adaptive_rejection_sampler_template(object):
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
        self.D =  np.array([self.X[0],self.X[self.b/2],X[self.b-1]])
        self.loga = np.sum(np.log(X))

    def add_to_D(self):
        out = np.array([],dtype=float64)
        new = np.random.gamma(10)
        counter = 0
        for i in self.D:
            if i <= new:
                out[counter] = i
            else:
                out[counter] = new
            counter = counter + 1
            
X = np.random.gamma(10,size=11)
lower_bound = 0.
upper_bound = np.inf
instance= adaptive_rejection_sampler_template(lower_bound,upper_bound,X)
print instance.D
instance.add_to_D()

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