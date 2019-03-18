### Utility Functions
from sklearn.preprocessing import normalize
import numpy as np
# Converts Distances to Kernels
def dist_to_kernel(D):
    
    gamma = 1.0 / np.amax(D)
    S = np.exp(-D * gamma)
    
    return S

# Computes Gram Matrix from Kernel
def GramMatrix(K):
    N_one = np.ones(K.shape) / K.shape[0] 
    K_bar = K - np.dot(N_one, K) - np.dot(K, N_one) + np.dot(np.dot(N_one, K), N_one)
    
    return K_bar


def GeneFilter(X, pct_dropout_min, pct_dropout_max):
    
    print("Applying Gene Filter . . .")

    dropouts = (np.sum(X == 0, axis = 1, keepdims = False)/X.shape[1]) * 100
    gene_filter = (dropouts < pct_dropout_max) & (dropouts > pct_dropout_min)

    print("Shape of Expression Matrix after Gene Filtering: ", X[gene_filter, :].shape)

    return X[gene_filter, :]


def FeatureNormalization(X, norm):
# X shape is (features, samples)
# norm is either l2, mean
# Returns a matrix X_nrm containing normalized values
    print ("Applying " + norm + " Normalization . . .")
    
    if (norm == "l2"):
        X_nrm = normalize(X.T, axis=0)
        X_nrm = X_nrm.T
        
    elif (norm == "mean"):
        mu = np.mean(X, axis=1, keepdims=True)
        sd = np.std(X, axis=1, keepdims=True)
        X_nrm = (X - mu)/sd
        
    elif (norm == "norm6"):
        min_f = np.amin(X, axis=1, keepdims=True)
        y = np.log(X + np.absolute(min_f) + 1)
        max_all = np.amax(y)
        X_nrm = y/max_all
        
    return X_nrm