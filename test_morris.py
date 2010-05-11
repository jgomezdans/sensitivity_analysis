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
    B_star = generate_trajectories ( numpy.array([1./3, 1./3]), \
                k, delta )
    pdb.set_trace()

def test_morris ():
    # A test of the Morris SA scheme
    # Use Sobol's g function, with 6 parameters
    # Sensitivity of each parameter is
    # inverseley proportional to parameter value
    a = numpy.array([78., 12., 0.5, 2., 97, 33])
    p = 4
    delta = 2/3.
    k = 2
    # drange defines the staring points for all
    # trajectories
    drange = numpy.arange ( 0, 4./3, 1./3)
    B_star = []
    # Create all trajectories. Define starting point
    # And calculate trajectory
    for i in itertools.product( drange, drange, drange, \
                                drange, drange, drange ):
        B_star.append (generate_trajectories ( numpy.array(i), \
            k, delta ) )
    # B_star contains all our trajectories
    B_star = numpy.array ( B_star )
    # Next stage: carry out the sensitivity analysis

    ee = [ k*([],)]
    for i in xrange(B_star.shape[0]):
        #for each trajectory, calculate the value of the model
        # at the starting point
        x0 = B_star[i,0,:]
        g_pre = sobol ( x0, a )
        for j in xrange(1, 7):
            # Get the new point in the trajectory
            x = B_star[i, j, :]
            #... and calculate the model output
            g = sobol( x, a )
            #store the difference. There's a denominator term here
            ee[numpy.nonzero(B_star[i,j,:]-B_star[i,j-1,:])[0]].append( g-g_pre )
            # Store the current value as the previous for the next
            # displacement along the trajectory
            g_pre = g
    # ee contains the distribution. Means and so on
    ee = numpy.array(ee)
    pdb.set_trace()

if __name__=="__main__":
    test_generate_trajectories()
    test_morris()
