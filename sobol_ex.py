# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from morris import *

def sobol ( x, a ):
    """
    The Sobol :math:`$g$` function  from Sobol 1990. Requires two parameters,
    """
    g = (np.abs(4.0*x -2 ) + a)/(1+a)
    return g.prod()

a = np.array([78., 12., 0.5, 2., 97, 33])
p = 4
delta = 2/3.
k = 6
# drange defines the staring points for all
# trajectories
num_traj = 100
drange = np.arange ( 0, 4./3, 1./3)
fsobol = lambda x, *a: sobol (x, np.array(a))
for t in [1]:#xrange(100):
    #(mu_star, mu, sigma) = sensitivity_analysis ( p, k, delta, \
                        #num_traj, drange, r=4, \
                        #func=fsobol, args=(a), sampling="campolongo" )

    (mu_star, mu, sigma) = sensitivity_analysis ( p, k, num_traj, \
        r=4, sampling="campolongo", \
        func=fsobol, args=(a) )

    #Labels=[ r'$x_{%1d}$'%i for i in xrange(1,7)]

    #plt.bar(np.arange(6)+.5,mu_star, width=0.5, fc='0.8', \
                        #label=r'$\mu_{i}^{*}$')
    #plt.ylabel(r'$\mu^{*}$')
    #plt.legend(loc='upper left', fancybox=True, shadow=True )
    #plt.twinx()
    #plt.bar(np.arange(6)+1,1./a, width=0.5, fc='0.4', \
                        #label=r'$1/a_{i}$')
    #plt.legend(loc='upper right', fancybox=True, shadow=True )
    #plt.ylabel(r'$1/a_{i}$')
    #ax = plt.gca()
    #ax.xaxis.set_ticks ( range(1,7))
    #ax.xaxis.set_ticklabels ( Labels )
    #plt.xlabel(r'Parameter')
    #plt.title('Example of sensitivity analysis of Sobol $g$ function.')
    #plt.show()

    mu_star = np.array ( mu_star )
    print mu_star.argsort()
    print (1./a).argsort()