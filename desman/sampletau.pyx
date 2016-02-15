"""
sampletau.pyx

"""

import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern int c_sample_tau (long *anTau, double* adPi, double *adEta, long* anVariants, int nV, int nG, int nS)

@cython.boundscheck(False)
@cython.wraparound(False)
def sample_tau(np.ndarray[long, ndim=3, mode="c"] tau not None, np.ndarray[double, ndim=2, mode="c"] pi not None, 
               np.ndarray[double, ndim=2, mode="c"] eta not None,
               np.ndarray[long, ndim=3, mode="c"] variants not None):
    """
    sample_tau (tau, pi, eta)
    Takes numpy arrays tau,pi, eta as input, samples tau one position and genome at a time
    param: tau -- a 3-d numpy array of np.int VXGX4
    param: pi -- strain frequencies SXG
    param: eta -- error rates 4X4
    param: variants - variant frequencies VXSX4
    """
    cdef int nV, nG

    nV = tau.shape[0]
    nG = tau.shape[1]
    nS = pi.shape[0]
    
    nchange = c_sample_tau (&tau[0,0,0], &pi[0,0], &eta[0,0], &variants[0,0,0],nV, nG, nS)

    return nchange