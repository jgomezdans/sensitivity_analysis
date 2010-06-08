# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from morris import *
import pdb

def sobol ( x, a ):
    """
    The Sobol :math:`$g$` function  from Sobol 1990. Requires two parameters,
    """
    g = (numpy.abs(4.0*x -2 ) + a)/(1+a)
    return g.prod()

a = np.array([78., 12., 0.5, 2., 97, 33])
p = 4
delta = 2/3.
k = 6
# drange defines the staring points for all
# trajectories
num_traj = 50
drange = numpy.arange ( 0, 4./3, 1./3)
fsobol = lambda x, *a: sobol (x, numpy.array(a))

for p_levels in xrange( 2, 12, 2):
    print p_levels
    for r_levels in xrange( 3, 7 ):
        print "\t", r_levels, 
        p = p_levels
        delta = p/(2*(p-1))
        k = 6
        num_traj = 100
        (mu_star, mu, sigma) = sensitivity_analysis ( p, k, delta, \
                        num_traj, drange, r=r_levels, \
                        func=fsobol, args=(a), sampling="campolongo" )
        mu_star = np.array ( mu_star )
        print numpy.argsort(-mu_star).argsort()
        try:
            ranker = np.vstack ([ ranker, ((-1*mu_star).argsort().argsort() + 1)])
        except NameError:
            ranker = ((-1*mu_star).argsort().argsort() + 1).T

