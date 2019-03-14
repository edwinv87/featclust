import numpy as np

from PairwiseDists import euclidean, pearsons, spearmans, linear_kernel

from UtilFunctions import dist_to_kernel, GramMatrix


def PCA(X, dist_or_kernel = 'linear', n_comp = 1):

    if (dist_or_kernel == 'linear'):
        K = GramMatrix(linear_kernel(X.T))

    elif (dist_or_kernel == 'spearmans'):
        K = GramMatrix(dist_to_kernel(spearmans(X.T)))

    elif (dist_or_kernel == 'euclidean'):
        K = GramMatrix(dist_to_kernel(euclidean(X.T)))

    elif (dist_or_kernel == 'pearsons'):
        K = GramMatrix(dist_to_kernel(pearsons(X.T)))

    E_vec, E_val, _ = np.linalg.svd(K)

    #print(E_val)
    #print(E_vec[:,0])
    # Sort Eigenvector 
    idx = np.argsort(E_val, kind = 'mergesort')
    idx = np.flip(idx)
    E_val = E_val[idx]
    E_vec = E_vec[:, idx]

    # Remove zero eigenvectors
    idx2 = E_val > 0
    E_val = E_val[idx2]
    E_vec = E_vec[:, idx2]
    print("Maximum components possible = ", E_val.size)

    # Scale eigenvectors so that np.dot(D[:,0].T, D[:, 0]) * E[0] = 1

    E_val = np.reshape(E_val, [1, E_val.size])
    # E_vec = E_vec / np.linalg.norm(E_vec, axis = 0, keepdims = True)
    E_vec = E_vec / np.sqrt(E_val)
    X_red = np.dot(E_vec[:, 0:n_comp].T, K)

    # print(E_vec[:,0])
    # print(X_red)
    return X_red
