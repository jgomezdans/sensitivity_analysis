# -*- coding: utf-8 -*-
import numpy

def sobol ( x, a ):
    """
    The Sobol :math:`$g$` function  from Sobol 1990. Requires two parameters, 
    """
    g = (numpy.abs(4.0*x -2 ) + a)/(1+a)
    return g.product()

def generate_trajectories ( x0, r, p, delta ):
    """
    Generate Morris trajectories to sample parameter space
    """
    k = x0.shape[0]

    if p%2 != 0:
        raise ValueError, "p number has to be even!"
    signo = numpy.random( k )
    signo = numpy.where (signo>0.5, 1, -1)
    D = numpy.matrix ( numpy.diag ( signo ) )
    D = numpy.matrix([1,0,0, -1]).reshape((k,k))
    P = numpy.zeros((k,k))
    pr = numpy.random.permutation ( k )
    for i in xrange(k):
        P[i, pr[i]] = 1
    P = numpy.matrix( P )
    B = numpy.matrix(numpy.tri(k+1, k, k=-1))
    J = numpy.matrix(numpy.ones ( (k+1, k)))
    B_star = (((2.0*B - J)*D + J)*(delta/2.) + J*x0)*P
    return B_star

test generate_trajectories ():
    k = 2
    p = 4
    delta = 2/3.
    generate_trajectories ( numpy.array([1./3, 1./3]) )
