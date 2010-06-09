# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from morris import *
import pdb

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
num_traj = 50
TDCC = np.zeros ( (3, 7) )#plevels, r_levels
P_VALS = np.zeros ( (3, 7) )#plevels, r_levels
fsobol = lambda x, *a: sobol (x, np.array(a))
i1 = 0
j1 = 0
for p_levels in xrange( 4, 10, 2):
    print p_levels
    for r_levels in xrange( 3, 10 ):
        j1 = 0
        for n_tries in xrange ( 50 ):
            print "\t", r_levels,
            (mu_star, mu, sigma) = sensitivity_analysis ( p, k, num_traj, \
                r=r_levels, sampling="campolongo", \
                func=fsobol, args=(a) )
            mu_star = np.array ( mu_star )
            try:
                ranker = np.vstack ([ ranker, \
                        ((-1*mu_star).argsort().argsort() + 1)])
            except NameError:
                ranker = ((-1*mu_star).argsort().argsort() + 1).T
        (tdcc, p_value ) = top_down_concordance ( ranker.T )
        del ranker
        TDCC[ i1, j1 ] = tdcc
        P_VALS[ i1, j1 ] = p_value
        j1 += 1
    i1 += 1