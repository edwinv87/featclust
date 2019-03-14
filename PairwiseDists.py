import numpy as np

from sklearn.metrics import pairwise_distances
from scipy.stats import spearmanr

### Pairwise Distances 

def spearmans(X):
    """
    X is a (n x d) matrix where rows are samples
    Returns D which is a (n x n) matrix of distances between samples
    
    """
    D, _ = spearmanr(X, axis = 1)
    D = 1 - D
    
    return D


def pearsons(X):
    """
    X is a (n x d) matrix where rows are samples
    Returns D which is a (n x n) matrix of distances between samples
    
    """
    D = pairwise_distances(X, metric = "correlation")
    
    return D
    
    
def euclidean(X):
    """
    X is a (n x d) matrix where rows are samples
    Returns D which is a (n x n) matrix of distances between samples
    
    """
    D = pairwise_distances(X, metric = "euclidean", n_jobs = -1)
    
    return D



### Kernels

def linear_kernel(X):
    """
    X is a (n x d) matrix where rows are samples
    Returns K which is a (n x n) kernel matrix
    
    """    
    K = np.dot(X, X.T)
    
    return K



