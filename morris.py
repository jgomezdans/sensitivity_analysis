# -*- coding: utf-8 -*-
import numpy as np
import itertools

#import pdb

def product(*args, **kwds):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = range(r)
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)


def generate_trajectory ( x0, p, delta ):
    """
    Generate Morris trajectories to sample parameter space

    :param x0: Initial trajectory location
    :param p: Number of quantisation levels for parameter space
    :param delta: The delta parameter from Saltelli et al.
    """
    k = x0.shape[0]

    if p % 2 != 0:
        raise ValueError, "p number has to be even!"
    signo = np.random.rand( k )
    signo = np.where (signo>0.5, 1, -1)
    D = np.matrix ( np.diag ( signo ) )
    #D = np.matrix([1,0,0, -1]).reshape((k,k))
    P = np.zeros((k, k))
    pr = np.random.permutation ( k )
    for i in xrange(k):
        P[i, pr[i]] = 1
    P = np.matrix( P )
    B = np.matrix(np.tri(k+1, k, k=-1))
    J = np.ones ( (k+1, k))
    B_star = (((2.0*B - J)*D + J)*(delta/2.) + J*x0)*P
    return B_star

######def campolongo_sampling ( b_star, r ):
    ######"""
    ######The campolongo sampling strategy, a brute-force search to find
    ######a set of r trajectories that would enable the best possible
    ######sampling of parameter space.

    ######My implementation is impractical as of yet!

    ######@param b_star: a (num_traj, k+1, k) trajectory matrix of elemental effects. A set of r that maximise parameter space exploration will beh chosen.
    ######"""
    #######import math
    ######import math
    ######num_traj = b_star.shape[0]
    ######k = b_star.shape[2]


    ####### Precalculate distances between all pairs of trajectories
    ######traj_distance = {}
    ######for ( m, l ) in product(range(num_traj), range(num_traj)):
        ######for ( i, j ) in product ( range(k), range(k) ):
            ######A = [ (b_star[m, i, z] - b_star[l, j, z])**2 \
                                    ######for z in xrange(k) ]
            ####### A will always be >0, so no need for sqrt
            ######traj_distance[ ( m, l ) ] = sum( A )#math.sqrt (sum(A))
            
        
    ####### Calculate aggregated distances by groups of trajectories
    ######selected_trajectories = list(([],)*8)
    ######for batches in xrange(8):
        ######traj_start = ( num_traj/8. )*batches
        ######traj_end = (num_traj/8.)*(batches+1)
        ######cnt = 0
        ######max_dist = 0.
        ######for h in combinations (range(traj_start, traj_end), r):
            ######cnt += 1
            ######accum = 0
            ######for (m,l) in combinations (h, 2):
                ######accum += traj_distance[ ( m, l ) ]
            ######if max_dist < accum:
                ######selected_trajectories[batches] =  h
                ######max_dist = accum
        
    ######selected_trajectories = np.array ( selected_trajectories ).flatten()
    ######cnt = 0
    ######traj = []
    ######distance = []
    ####### Now, we can pick and mix the trajectories from the best sets
    ######for h in combinations (selected_trajectories, r):
        ######cnt += 1
        ######accum = 0
        ######for (m,l) in combinations (h, 2):
            ######accum += traj_distance[ ( m, l ) ]
        ######if max_dist < accum:
            
            ######traj.append( h )
            ######distance.append ( accum )
            ######max_dist = accum

    ######distance = np.array ( distance )
    ######i = distance.argsort()
    ######s = np.unique ( np.array ( traj) [i][-(r+1):] )
    ######return b_star[ s, :, :]
            
def sensitivity_analysis ( p, k, delta, num_traj, drange, \
                           func, args=(), r=None, \
                           sampling="Morris" ):
    """
    Carry out a sensitivity analysis using the Morris approach.

    @param p: The :math:`p` parameter: parameter space quantisation.
    @param k: Number of parameters.
    @param num_traj: Number of trajectories to calculate (Morris method)
    @param drange: Ranges to be used.
    @param func: Model function
    @param args: extra arguments to model function
    @param r: Campolongo;s trajectories (:math:`r<num_traj`)
    @param sampling: Sampling type. Either "Morris" or "Campolongo"
    """
    if sampling != "Morris":
        if sampling.lower() != "campolongo":
            raise ValueError, "For Campolongo scheme, r >0"
        if r == 0:
            raise ValueError, "Need a subset of chains"
    B_star = []
    # Create all trajectories. Define starting point
    # And calculate trajectory
    counter = 0
    for i in product( drange, drange, drange, \
                                drange, drange, drange ):
        if np.random.rand() > 0.5:
            B_star.append (generate_trajectory ( np.array(i), \
                k, delta ) )
            counter += 1
            if counter > num_traj: break
    # B_star contains all our trajectories
    B_star = np.array ( B_star )
    # Next stage: carry out the sensitivity analysis
    if sampling != "Morris":
        B_star = campolongo_sampling ( B_star, r )
    ee = [ [] for i in xrange(k) ]
    for i in xrange(B_star.shape[0]):
        #for each trajectory, calculate the value of the model
        # at the starting point
        x0 = B_star[i, 0, :]
        g_pre = func ( x0, *args )
        for j in xrange(1, 7):
            # Get the new point in the trajectory
            x = B_star[i, j, :]
            #... and calculate the model output
            g = func( x, *args )
            #store the difference. There's a denominator term here
            idx = np.nonzero(B_star[i, j, :] - \
                    B_star[i, j-1, :])[0]
            ee[idx].append( g-g_pre )
            
            # Store the current value as the previous for the next
            # displacement along the trajectory
            g_pre = g
    # ee contains the distribution. Means and so on
    #pdb.set_trace()
    E = [ np.array(x) for x in ee]
    mu_star =[ np.abs(u).mean() for u in E]
    mu =[ u.mean() for u in E]
    sigma =[ u.std() for u in E]
    return ( mu_star, mu, sigma )

def top_down_concordance ( replicates_matrix ):
    """
    Top-down concordance coefficient (TDCC) from Iman and Conover.

    :parameter ranks: A matrix (num_parameters, replicate) of parameter
                      rankings
    """

    ( k, replicates ) = replicates_matrix.shape
    ranks = np.array([ replicates_matrix[:, i] for i in xrange( replicates )])
    nsa = replicates
    # Calculate savage scores... Groar!!!
    ss = np.zeros_like ( ranks )
    ss = [ np.sum(1./np.arange( ri, k+1)) for ri in ranks.flat ]
    numerator = (ss.sum(axis=0)**2).sum() -nsa*nsa*k
    denominator = nsa*nsa*( k-(1./np.arange(1, k+1)).sum() )
    
    tdcc = numerator / denominator
    p_value = nsa * ( k - 1 ) * tdcc
    return (tdcc, p_value )

def cliqdistances( cliq, dist ):
    return sorted( [dist[j,k] for j in cliq  for k in cliq if j < k], reverse=True )

def maxarray2( a, n ):
    """ -> max n [ (a[j,k], (j,k)) ...]  j <= k, a symmetric """
    jkflat = np.argsort( a, axis=None )[:-2*n:-1]
    jks = [np.unravel_index( jk, a.shape ) for jk in jkflat]
    return [(a[j,k], (j,k)) for j,k in jks if j <= k] [:n]

def _str( iter, fmt="%.2g" ):
    return " ".join( fmt % x  for x in iter )

#...............................................................................

def maxweightcliques( dist, nbest, r, verbose=10 ):

    def cliqwt( cliq, p ):
        return sum( dist[c,p] for c in cliq )  # << 0 if p in c

    def growcliqs( cliqs, nbest ):
        """ [(cliqweight, n-cliq) ...] -> nbest [(cliqweight, n+1 cliq) ...] """
            # heapq the nbest ? here just gen all N * |cliqs|, sort
        all = []
        dups = set()
        for w, c in cliqs:
            for p in xrange(N):
                    # fast gen [sorted c+p ...] with small sorted c ?
                cp = c + [p]
                cp.sort()
                tup = tuple(cp)
                if tup in dups:  continue
                dups.add( tup )
                all.append( (w + cliqwt(c, p), cp ))
        all.sort( reverse=True )
        if verbose:
            print "growcliqs: %s" % _str( w for w,c in all[:verbose] ) ,
            print " best: %s" % _str( cliqdistances( all[0][1], dist )[:10])
        return all[:nbest]

    dist.flat[::dist.shape[0]+1] = -1e10
    #np.fill_diagonal( dist, -1e10 )  # so cliqwt( c, p in c ) << 0
    C = (r+1) * [(0, None)]  # [(cliqweight, cliq-tuple) ...]
        # C[1] = [(0, (p,)) for p in xrange(N)]
    C[2] = [(w, list(pair)) for w, pair in maxarray2( dist, nbest[2] )]
    for j in range( 3, r+1 ):
        C[j] = growcliqs( C[j-1], nbest[j] )
    return C
    
def clique_percolation ( b_star, r ):
    from  __future__ import division
    N = 50
    r = 2  # max clique size
    nbest = 3
    verbose = 10
    seed = 1
    exec "\n".join( sys.argv[1:] )  # N= ...
    #np.random.seed(seed)
    nbest = [0, 0, N//2] + (r - 2) * [nbest]  # ?

    print "  N=%d  r=%d  nbest=%s"  % ( N, r, nbest)

        # random graphs w cluster parameters ?
    dist = np.random.exponential( 1, (N,N) )
    dist = (dist + dist.T) / 2
    for j in range( 0, N, r ):
        dist[j:j+r, j:j+r] += 2  # see if we get r in a row
    # dist = np.ones( (N,N) )

    cliqs = maxweightcliques( dist, nbest, r, verbose )[-1]  # [ (wt, cliq) ... ]

    print "Clique weight,  clique,  distances within clique"
    print 50 * "-"
    for w,c in cliqs:
        print "%5.3g  %s  %s" % (
            w, _str( c, fmt="%d" ), _str( cliqdistances( c, dist )[:10]))

    
def campolongo_sampling2 ( b_star, r ):
    """
    The campolongo sampling strategy, a brute-force search to find
    a set of r trajectories that would enable the best possible
    sampling of parameter space.

    My implementation is impractical as of yet!

    @param b_star: a (num_traj, k+1, k) trajectory matrix of elemental effects. A set of r that maximise parameter space exploration will beh chosen.
    """
    #import math
    import math
    from  __future__ import division
    num_traj = b_star.shape[0]
    k = b_star.shape[2]

    #---------------------------------------------------------
    # Precalculate distances between all pairs of trajectories
    #---------------------------------------------------------
    traj_distance = {}
    for ( m, l ) in product(range(num_traj), range(num_traj)):
        for ( i, j ) in product ( range(k), range(k) ):
            A = [ (b_star[m, i, z] - b_star[l, j, z])**2 \
                                    for z in xrange(k) ]
            # A will always be >0, so no need for sqrt
            traj_distance[ ( m, l ) ] = sum( A )#math.sqrt (sum(A))

    N = num_traj
    r = 8  # max clique size
    nbest = 10
    verbose = 10
    seed = 1
    nbest = [0, 0, N//2] + (r - 2) * [nbest]  # ?

    print "  N=%d  r=%d  nbest=%s"  % ( N, r, nbest)

        # random graphs w cluster parameters ?
    
    dist = (traj_dist + traj_dist.T) / 2
    for j in range( 0, N, r ):
        dist[j:j+r, j:j+r] += 2  # see if we get r in a row
    # dist = np.ones( (N,N) )

    cliqs = maxweightcliques( dist, nbest, r, verbose )[-1]  # [ (wt, cliq) ... ]

    print "Clique weight,  clique,  distances within clique"
    print 50 * "-"
    for w,c in cliqs:
        print "%5.3g  %s  %s" % (
            w, _str( c, fmt="%d" ), _str( cliqdistances( c, dist )[:10]))

    

