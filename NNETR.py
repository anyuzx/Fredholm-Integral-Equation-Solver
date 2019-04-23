import numpy as np
import scipy

# core algorithm of non-negative equality Tikhonov regularization (NNETR)
def NNETR(K, f, Delta, epsilon, alpha):
    # the first step
    A_nn = np.vstack((K, alpha * np.identity(K.shape[1])))
    b_nn = np.hstack((f, np.zeros(K.shape[1])))
    
    # the second step
    A_nne = np.vstack((epsilon * A_nn, np.full(A_nn.shape[1], 1.0)))
    b_nne = np.hstack((epsilon * b_nn, 1.0))

    # Use NNLS solver provided by scipy
    sol, residue = scipy.optimize.nnls(A_nne, b_nne)

    # solution should be divided by Delta (grid size)
    sol = sol/Delta
    return sol, residue
