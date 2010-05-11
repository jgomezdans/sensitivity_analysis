#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pdb
from morris import *

def sobol ( x, a ):
    """
    The Sobol :math:`$g$` function  from Sobol 1990. Requires two parameters,
    """
    g = (numpy.abs(4.0*x -2 ) + a)/(1+a)
    return g.prod()



def test_generate_trajectory ():
    k = 2
    p = 4
    delta = 2/3.
    B_star = generate_trajectory ( numpy.array([1./3, 1./3]), \
                k, delta )
    print B_star

def test_morris ():
   
    # A test of the Morris SA scheme
    # Use Sobol's g function, with 6 parameters
    # Sensitivity of each parameter is
    # inverseley proportional to parameter value
    a = numpy.array([78., 12., 0.5, 2., 97, 33])
    p = 4
    delta = 2/3.
    k = 6
    # drange defines the staring points for all
    # trajectories
    num_traj = 100
    drange = numpy.arange ( 0, 4./3, 1./3)
    fsobol = lambda x, *a: sobol (x, numpy.array(a))
    (mu_star, mu, sigma) = sensitivity_analysis ( p, k, delta, \
        num_traj, drange, \
        func=fsobol, args=(a), r=10, \
        sampling="Campolongo" )
    print mu_star
    print sigma
    print 1./a
    pdb.set_trace()
if __name__=="__main__":
    test_generate_trajectory()
    test_morris()
