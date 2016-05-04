# -*- coding: utf-8 -*-
"""
Created on Mon May  2 19:44:58 2016
The idea here is that we should do this a little more thoughtfully 
so that calculations are relatively easy.
@author: brianlibgober
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stat
from scipy.stats import gaussian_kde

from scipy.stats import kstest
from matplotlib.backends.backend_pdf import PdfPages

#%%
class Interval:
    def __init__(self,lower,upper):
        self.tuple = (lower,upper)
        self.lower = lower
        self.upper = upper
        self.length = upper-lower
        
    def interior_contains(self,x):
        return self.lower < x < self.upper
    
    def closure_contains(self,x):
        return self.lower <= x <= self.upper
        
    def __repr__(self):
        return "(" + str(self.lower) + "," + str(self.upper) + ")" 
        
    
class Segment(Interval):
    """
    Container for line segment information.
    @ lower : lower bound of segment
    @ upper : upper bound of segment
    @ point : location of point anywhere on the segment
    @ slope : slope of the line segment
    
    Methods include calculating the value
    """
    def __init__(self,lower,upper,point,slope):
        Interval.__init__(self,lower,upper)
        self.point = point
        self.slope = slope
        self.mass = None

    @property
    def y_intercept(self): 
        return  -1*self.slope*self.point[0] + self.point[1] 
        
    @property
    def pointslope(self):
        b = self.y_intercept  
        if b>0:
            return  str(self.slope)+"x+" + str(b)
        else:
            return str(self.slope)+"x+" + str(b)
         
    def __repr__(self):
        return  "(" + self.pointslope+  ")" \
            "*I{[" + str(self.lower) + "," + str(self.upper) + "]}"\
            
    def value(self,x):
        if not self.closure_contains(x):
            print "X not in closure of segment"
            return None
        else:
            return self.slope*(x-self.point[0]) + self.point[1]

    @property
    def uppervalue(self):
        return self.slope*(self.upper-self.point[0]) + self.point[1]
    
    @property
    def lowervalue(self):
        return self.slope*(self.lower-self.point[0]) + self.point[1]

    @property 
    def maxvalue(self):
        return max(self.lowervalue,self.uppervalue)

    
self = Segment(1,2,(0,0),2)
print self.value(2)
print self
print "Upper value", self.uppervalue
print "Lower value", self.lowervalue
print "Max value", self.maxvalue

# %% Define a Partition class


class Partition:
    """
    A container for a partition of the real line.
    We assume it is made up of half closed intervals like [a,b)
    
    Input:
    @ knots : locations in the interval
    @ xlb : lower bound
    @ xub : upper bound
    """


    def __init__(self,knots, xlb=-np.inf,xub=np.inf):
        self.xlb = xlb
        self.xub = xub
        self.domain = Interval(xlb,xub)
        self._knots = np.unique(np.concatenate(([xlb],np.sort(knots),[xub])))
        self._intervals = np.array(
            [Interval(self._knots[i],self._knots[i+1])
                    for i in xrange(self._knots.size-1)])
    
    @property
    def knotcount(self):
        return self._knots.size
        
    @property
    def knots(self):
        return self._knots
    
    @property
    def interior_knots(self):
        return self._knots[1:-1]
        
    def __repr__(self):
        return self._intervals.__repr__()
        
    def insert_knot(self,x):
        """
        Takes a number and if it not already a knot in the partition adds it.
        Also updates the intervals
        
        Returns the zero indexed location of new knot in the partition. 
        """
        if (x not in self.knots) and self.domain.closure_contains(x):
            interval_number = 0
            for interval in self._intervals:
                if interval.interior_contains(x):
                    #replace the old interval with one broken around x
                    #add another from x to old upper
                    self._knots = np.insert(self._knots,interval_number+1,x)
                    self._intervals[interval_number] = Interval(interval.lower,x)
                    self._intervals = np.insert(self._intervals, \
                       interval_number+1,
                       Interval(x,interval.upper))
                    break
                else:
                    interval_number += 1
            return interval_number + 1
        else:
            return None
        
self = Partition([1,2])
self.knotcount
self.interior_knots
print self
self.insert_knot(5)
print self
self.insert_knot(10)
# %%
class Hull():
    """
    A container for a hull class.  This is like a collection of segments.
    
    Input:\n    
    @ knots : segment edgepoints, excluding the  boundaries
    @ values: values of the hull at each input knot
    @ slopes  : slopes for the initialized segments
    @ xlb : lower bound
    @ xub : upper bound
    """
    def __init__(self,knots,values, slopes, xlb=-np.inf,xub=np.inf):
        self.xlb = xlb
        self.xub = xub
        self.domain = Interval(xlb,xub)
        self._values = np.concatenate(([xlb],values,[xub]))
        self._knots = np.unique(np.concatenate(([xlb],np.sort(knots),[xub])))
        self._segments = np.array(Segment(lower=self._knots[0],
                                          upper =self._knots[1],
                                    point=(self._knots[1],values[0]),
                                    slope = slopes[0])
                                    )
        for i in xrange(1,self._knots.size-2):
            self._segments = np.append(self._segments,
                Segment(lower=self._knots[i],
                        upper=self._knots[i+1],
                        point =(self._knots[i],values[i-1]),
                        slope=slopes[i])
                        )
        self._segments = np.append(self._segments,
                Segment(lower=self._knots[self._knots.size-2],
                        upper=self._knots[self._knots.size-1],
                        point=( #want second to last _knot, -1
                        #and zero-indexing => -1 again
                        self._knots[self._knots.size-2],
                        #self._knots has two extra points over values
                        #and then there's zero indexing
                        values[self._knots.size-3]),
                        #there are one fewer slopes
                        #zero indexing implies we get this number +1th
                        slope=slopes[self._knots.size-2])
                )
                
    def __repr__(self):
        return self._segments.__repr__()
    
    def value(self,x):
        for segment in self._segments:
            if segment.closure_contains(x):
                return segment.value(x)
                
    def __getitem__(self,i):
        return self._segments[i]
    
    @property
    def interval_count(self):
        return self._knots.size-1
    
    @property
    def max_value(self):
        max_value = self._segments[0].maxvalue
        for i in xrange(1,self.interval_count):
            max_value = max(max_value,self._segments[1].maxvalue)
            if self._segments[i].slope <= 0:
                return max_value
            else:
                continue
        
                
            
self = Hull(knots=[1,2],values=[3,5],slopes=[1,0,-1],xlb=-5)     
print self    
print "Value at 0",self.value(0)
print "Value at 1",self.value(1)
print "Value at 2",self.value(2)


x  = np.arange(0,5,0.25)
plt.close("all")
plt.plot(x,np.vectorize(self.value)(x))
#%%
    
class ArsSampler():
    """
    This is an object that handles ARS sampling,
    It contains a partition and a hull.
    It also has methods for updating these and plotting them.
    
    @
    """
    def __init__(self,initial_knots,h,hprime,xlb=-np.inf,xub=np.inf,max_knots=50):
        self._X = Partition(initial_knots,xlb,xub)
        self._X = Partition(initial_knots,xlb,xub)
        self.h  = h
        self.hprime = hprime
        self._X.hknots = h(self._X.interior_knots) #hknots is defined on the interior knots
        self._X.hprimeknots =  hprime(self._X.interior_knots) #hprimeknots is not
        self.max_knots = max_knots
        self.Z = np.array([])
        self.hZ = np.array([])
        for i in xrange(self._X.interior_knots.size-1):
            xi = self._X.interior_knots[i]
            xip1 = self._X.interior_knots[i+1]
            hxi = self._X.hknots[i]
            hxip1 = self._X.hknots[i+1]
            hxi_prime = self._X.hprimeknots[i]
            hxip1_prime = self._X.hprimeknots[i+1]
            zi = (hxip1 - hxi - xip1*hxip1_prime + xi*hxi_prime)/(
            hxi_prime - hxip1_prime)
            self.Z = np.append(self.Z,zi)
            self.hZ = np.append(self.hZ,hxi_prime*(zi-xi)+hxi)         
        self._UpperHull = Hull(knots=self.Z,values=self.hZ,
                               slopes=self._X.hprimeknots,
                               xlb=xlb,
                               xub = xub)  
            

    @property
    def interval_count(self):    
        return
        
    def z(self,xi,hxi,hxi_prime,xip1,hxip1,hxip1_prime):
        return xi+(hxi - hxip1 + hxip1_prime*(xip1  - xi))/(
        hxip1_prime - hxi_prime)
    
    
    def show(self):
        import matplotlib as plot
        self._X.interior_knots[-1] - self._X.interior_knots[0]
        xpoints  = np.arange(1.1*self._X.interior_knots[0],
                       1.1*self._X.interior_knots[-1],
                        min(
                        [(self._X.interior_knots[-1] - self._X.interior_knots[0])/100,0.25]
                        ))
        np.vectorize(self._UpperHull.value)
        plt.plot(xpoints,self.h(xpoints),"b")
        plt.plot(self._X.interior_knots,self._X.hknots,"rs")
        hullfunc = np.vectorize(self._UpperHull.value)
        plt.plot(self._UpperHull._knots,hullfunc(self._UpperHull._knots),"g^")
        plt.plot(xpoints,np.vectorize(self._UpperHull.value)(xpoints),"black")    
    
    
    def insert_value(self,newvalue,evaluated=None):
        #try to insert into partition, if already there or outside domain
        #then will return none
        tmp = self._X.insert_knot(newvalue)
        if tmp is None:
            return None
        interior_index = tmp-1
        if evaluated is None:
            evaluated = self.h(newvalue)
        self._X.hknots = np.insert(self._X.hknots,interior_index,evaluated)
        self._X.hprimeknots = np.insert(self._X.hprimeknots,interior_index,self.hprime(newvalue))
        #assuming that needs to be a point to the right,add it
        if 0 < interior_index < self._X.interior_knots.size-1:
            xim1 = self._X.interior_knots[interior_index-1]
            xi = self._X.interior_knots[interior_index]
            xip1 = self._X.interior_knots[interior_index+1]
            hxim1 = self._X.hknots[interior_index-1]
            hxi = self._X.hknots[interior_index]
            hxip1 = self._X.hknots[interior_index+1]
            hxim1_prime = self._X.hprimeknots[interior_index-1]
            hxi_prime = self._X.hprimeknots[interior_index]
            hxip1_prime = self._X.hprimeknots[interior_index+1]
            ziabove = self.z(xi,hxi,hxi_prime,xip1,hxip1,hxip1_prime)
            zibelow = self.z(xim1,hxim1,hxim1_prime,xi,hxi,hxi_prime)
            oldsegment_below = self._UpperHull._segments[interior_index-1]
            oldsegment_above= self._UpperHull._segments[interior_index]
            #reach into Upper Hull type
            self._UpperHull._knots = np.concatenate((self._UpperHull._knots[:interior_index],
                                                    [zibelow,ziabove],
                                                    self._UpperHull._knots[interior_index+1:]))

            self._UpperHull._segments = np.concatenate( (
                            self._UpperHull._segments[:interior_index-1],
                            [Segment(lower=oldsegment_below.lower,
                                     upper= zibelow,
                                     point = (xim1,hxim1),
                                    slope = oldsegment_below.slope),
                            Segment(lower=zibelow,
                                     upper=ziabove,
                                     point =(xi,hxi),
                                    slope = hxi_prime),
                            Segment(lower=ziabove,
                                     upper=oldsegment_above.upper,
                                     point = (xip1,hxip1),
                                    slope = oldsegment_above.slope)],
                             self._UpperHull._segments[interior_index+1:]
                            ))    

                            
        elif interior_index == 0:
            xi = self._X.interior_knots[0]
            xip1 = self._X.interior_knots[1]
            hxi = self._X.hknots[0]
            hxip1 = self._X.hknots[1]
            hxi_prime = self._X.hprimeknots[0]
            hxip1_prime = self._X.hprimeknots[1]
            ziabove = self.z(xi,hxi,hxi_prime,xip1,hxip1,hxip1_prime)
            oldsegment_above= self._UpperHull._segments[0]
            self._UpperHull._knots = np.concatenate(([self._UpperHull.xlb,
                                                      ziabove],
                                                    self._UpperHull._knots[1:]))
                                                    
            self._UpperHull._segments = np.concatenate( (
                            [Segment(lower=self._UpperHull.xlb,
                                     upper= ziabove,
                                     point = (xi,hxi),
                                     slope = hxi_prime),
                            Segment(lower=ziabove,
                                     upper=oldsegment_above.upper,
                                     point =oldsegment_above.point,
                                    slope = oldsegment_above.slope)],
                             self._UpperHull._segments[1:]
                            ))  


        else:
            end = self._X.interior_knots.size-1
            xim1 = self._X.interior_knots[end-1]
            xi = self._X.interior_knots[end]
            hxim1 = self._X.hknots[end-1]
            hxi = self._X.hknots[end]
            hxim1_prime = self._X.hprimeknots[end-1]
            hxi_prime = self._X.hprimeknots[end]
            zibelow = self.z(xim1,hxim1,hxim1_prime,xi,hxi,hxi_prime)
            oldsegment_below = self._UpperHull._segments[end-1]
            #reach into Upper Hull type
            self._UpperHull._knots = np.concatenate((self._UpperHull._knots[:end],
                                                    [zibelow,self._UpperHull.xub])
                                                    )
                                                    
            self._UpperHull._segments = np.concatenate( (
                            self._UpperHull._segments[:end-1],
                            [Segment(lower=oldsegment_below.lower,
                                     upper= zibelow,
                                     point = (xim1,hxim1),
                                    slope = oldsegment_below.slope),
                            Segment(lower=zibelow,
                                     upper=self._UpperHull.xub,
                                     point =(xi,hxi),
                                    slope = hxi_prime)]
                                    ) )


    @property
    def _masses(self):
        """
        If c = integral of exp(h), where h is the upper hull
        then K = c/exp(h*) where  h* is the max of the upper hull.
        Note that np.exp(1000)=inf, so if we want numerical stabilitiy
        We should divide everything through by this large value.
        """
        masses = np.array([])
        hullmax =  self._UpperHull.max_value
        for segment in self._UpperHull._segments:
            if segment.slope !=0:
                masses =  np.append(masses,(1./segment.slope)*\
                (np.exp(segment.uppervalue - hullmax) - \
                np.exp(segment.lowervalue - hullmax)))
            else:
                masses = np.append(masses,
                segment.length*np.exp(segment.lowervalue-hullmax))
        return masses
     
    @property
    def max_value(self):     
        return self._UpperHull.max_value        
        
    def _sample_hull(self):
        """
        Draws samples from the distribution whose log is given by the hull
        """
        U = np.random.random(1)[0]
        masses = self._masses
        K= masses.sum()
        cumulated = np.insert(np.cumsum(masses/K),0,0)
        for i in xrange(self._UpperHull.interval_count):
            if U < cumulated[i+1]:
                segment = self._UpperHull[i]
                m = segment.slope
                b = segment.y_intercept
                hstar = self.max_value
                if segment.lower != -np.inf:
                    if segment.slope != 0:
                        return segment.lower + 1./m*(
                        np.log(
                        1+(m*(U-cumulated[i])*K)/
                        np.exp(segment.lowervalue-hstar))
                            )
                    else:
                        return (
                        (U-cumulated[i])*K*np.exp(hstar-b)+segment.lower
                        )
                else:
                   return (np.log(U*m*K) -
                   (b-hstar))/m
        return "Fail"
        
    def _sample(self):
        """
        Draw a single sample from the distribution
        """
        x = self._sample_hull()
        hx = self.h(x)
        u = np.random.random()
        ACCEPTED = np.log(u) < hx - self._UpperHull.value(x)
        if ACCEPTED:
            return x,hx,ACCEPTED
        else:
            return x,hx,ACCEPTED
    
    def sample(self,size):
         samples = np.array([])
         while samples.size < size:
             x, hx, ACCEPTED = self._sample()
             if ACCEPTED:
                 samples = np.append(samples,x)
             if self._UpperHull._knots.size < self.max_knots + 2:
                 self.insert_value(newvalue=x,evaluated=hx)
         
         return samples
        
        
        
#%% #reinitialize empty class for prototyping, now we're inside the class
np.random.seed(1)
initial_knots = [-1,0,1]
h = np.vectorize(lambda x: (-1.*(x)**2)/2) #log of  a normal(0,1) pdf, up to proportionality constant
hprime = np.vectorize(lambda x: -x)
xlb = -np.inf
xub = np.inf
self = ArsSampler(initial_knots,h,hprime,xlb,xub)
samples = self.sample(10000)
xs = np.linspace(-5,5,200)
distro = stat.norm()
fig1 = plt.figure()
ax1 = fig1.add_subplot(2,2,1)
ax1.plot(xs,gaussian_kde(samples)(xs))
ax1.plot(xs,distro.pdf(xs))
kstest(samples,"norm")
#%%
#np.random.seed(1)
#for i in np.random.random(size=100):
#    self.insert_value(30*(i-0.5))
#self.show()
#func = np.vectorize(self._UpperHull.value)
#%%
#np.random.seed(2)
#samples = self.sample(10000)
#print samples
#plt.hist(samples)
#from scipy.stats import kstest
#kstest(samples,"norm")
#%% Gamma k, theta
k = 9; theta=0.5
distro = stat.gamma(a=k,scale=0.5)
h = np.vectorize(lambda x: (k-1)*np.log(x) -x/theta)
hprime = np.vectorize(lambda x: (k-1)/x -1/theta)
initial_knots = [1,4.5,7]
xlb = -np.inf
xub = np.inf
self = ArsSampler(initial_knots,h,hprime,xlb,xub)
np.random.seed(1)
samples = self.sample(10000)
xs = np.linspace(0,10,200)
ax2 = fig1.add_subplot(2,2,2)
ax2.plot(xs,gaussian_kde(samples)(xs))
ax2.plot(xs,distro.pdf(xs))
kstest(samples,distro.cdf)
#%% Truncated Normal, cut above -1
a=-1
b=np.inf
def h(x):
    if a <= x < b:
        return (-1.*(x)**2)/2
    else:
        return -np.inf

def hprime(x):
    if a <= x < b:
        return -1.*(x)
    else:
        return -np.nan

initial_knots = [-0.5,1,3]
xlb = a
xub = np.inf
distro = stat.truncnorm(a=a,b=np.inf)
np.random.seed(1)
self = ArsSampler(initial_knots,np.vectorize(h),np.vectorize(hprime),xlb,xub)
samples = self.sample(10000)
from scipy.stats import gaussian_kde
xs = np.linspace(-2,10,200)
ax3 = fig1.add_subplot(2,2,3)
ax3.plot(xs,gaussian_kde(samples)(xs))
ax3.plot(xs,distro.pdf(xs))
kstest(samples,distro.cdf)

#%% Bayesian Inference on  our conjugate distribution
np.random.seed(1)
X = np.sort(np.random.gamma(10,size=750))
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
self = ArsSampler(initial_knots,np.vectorize(h),np.vectorize(hprime),0,xub)
draws = self.sample(1000)
xs = np.linspace(9,11,200)
ax4 = fig1.add_subplot(2,2,4)
ax4.plot(xs,gaussian_kde(draws)(xs))
