# -*- coding: utf-8 -*-
import numpy as np
#from future import __division__

def cliqdistances( cliq, dist ):
    return sorted( [dist[j,k] for j in cliq  for k in cliq if j < k], reverse=True )

def maxarray2( a, n ):
    """ -> max n [ (a[j,k], (j,k)) ...]  j <= k, a symmetric """
    jkflat = np.argsort( a, axis=None )[:-2*n:-1]
    jks = [np.unravel_index( jk, a.shape ) for jk in jkflat]
    return [(a[j,k], (j,k)) for j,k in jks if j <= k] [:n]

def e_str( iter, fmt="%.2g" ):
    return " ".join( fmt % x  for x in iter )

#...............................................................................

def maxweightcliques( dist, nbest, r, N, verbose=10 ):

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
            print "growcliqs: %s" % e_str( w for w,c in all[:verbose] ) ,
            print " best: %s" % e_str( cliqdistances( all[0][1], dist )[:10])
        return all[:nbest]

    dist.flat[::dist.shape[0]+1] = -1e10
    #np.fill_diagonal( dist, -1e10 )  # so cliqwt( c, p in c ) << 0
    C = (r+1) * [(0, None)]  # [(cliqweight, cliq-tuple) ...]
        # C[1] = [(0, (p,)) for p in xrange(N)]
    C[2] = [(w, list(pair)) for w, pair in maxarray2( dist, nbest[2] )]
    for j in range( 3, r+1 ):
        C[j] = growcliqs( C[j-1], nbest[j] )
    return C