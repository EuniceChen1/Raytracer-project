
"""
A Genpolar task about generators

Eunice Chen 1/3/2014
"""
import numpy as np


def rtpairs(R,N):
    '''generates r,n pairs
        r and n are lists of the same amount of numbers 
    '''
    for i in range(len(R)):
	r=R[i]
	n=N[i]
        theta = 0.0 #starting with theta = 0
        for j in range(N[i]): 
            theta +=2*(np.pi)/N[i]
            yield R[i],theta


def rtuniform(n, rmax, m):
    '''generate r,n with uniform spacing
    '''
    R=np.arange(0,rmax,rmax/n)
    N=np.arange(1, n*m, m)
    return rtpairs(R, N)




