# -*- coding: utf-8 -*-
import numpy
import itertools
import pdb

def generate_trajectory ( x0, p, delta ):
    """
    Generate Morris trajectories to sample parameter space
    """
    k = x0.shape[0]

    if p%2 != 0:
        raise ValueError, "p number has to be even!"
    signo = numpy.random.rand( k )
    signo = numpy.where (signo>0.5, 1, -1)
    D = numpy.matrix ( numpy.diag ( signo ) )
    #D = numpy.matrix([1,0,0, -1]).reshape((k,k))
    P = numpy.zeros((k,k))
    pr = numpy.random.permutation ( k )
    for i in xrange(k):
        P[i, pr[i]] = 1
    P = numpy.matrix( P )
    B = numpy.matrix(numpy.tri(k+1, k, k=-1))
    J = numpy.ones ( (k+1, k))
    B_star = (((2.0*B - J)*D + J)*(delta/2.) + J*x0)*P
    return B_star

def campolongo_sampling ( b_star, r ):
    num_traj = b_star.shape[0]
    max_dist = 0.
    for h in itertools.combinations (range(num_traj), 4):
        for (m,l) in itertools.combinations (h, 2):
            accum = 0.
            for ( i, j ) in itertools.izip ( range(k), range(k) ):
                A = [ (b_star[m, i, z] - b_star[l, j, z])**2 \
                        for z in xrange(k) ]
                A = numpy.array( numpy.sqrt (A) ).sum()
                accum += A
        if max_dist < accum:
            selected_trajectories = h
    return b_star[ selected_trajectories, :, :]
            
def sensitivity_analysis ( p, k, delta, num_traj, drange, \
                           func, args=(), r=None, \
                           sampling="Morris" ):
    if sampling != "Morris":
        assert r != None
        raise ValueError, "For Campolongo scheme, r >0"
    B_star = []
    # Create all trajectories. Define starting point
    # And calculate trajectory
    for i in itertools.product( drange, drange, drange, \
                                drange, drange, drange ):
        B_star.append (generate_trajectory ( numpy.array(i), \
            k, delta ) )
    # B_star contains all our trajectories
    B_star = numpy.array ( B_star )
    # Next stage: carry out the sensitivity analysis
    if sampling != "Morris":
        B_star = campolongo_sampling ( B_star, r )
    ee = [ [] for i in xrange(k) ]
    for i in xrange(B_star.shape[0]):
        #for each trajectory, calculate the value of the model
        # at the starting point
        x0 = B_star[i,0,:]
        g_pre = func ( x0, *args )
        for j in xrange(1, 7):
            # Get the new point in the trajectory
            x = B_star[i, j, :]
            #... and calculate the model output
            g = func( x, *args )
            #store the difference. There's a denominator term here
            idx = numpy.nonzero(B_star[i, j, :] - \
                    B_star[i, j-1, :])[0]
            ee[idx].append( g-g_pre )
            
            # Store the current value as the previous for the next
            # displacement along the trajectory
            g_pre = g
    # ee contains the distribution. Means and so on
    pdb.set_trace()
    E = [ numpy.array(x) for x in ee]
    mu_star =[ numpy.abs(u).mean() for u in E]
    mu =[ u.mean() for u in E]
    sigma =[ u.std() for u in E]
    return ( mu_star, mu, sigma )

