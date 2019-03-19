# Import Packages
import pandas as pd
import numpy as np 

from sklearn.cluster import FeatureAgglomeration, AgglomerativeClustering
from sklearn.metrics import silhouette_score

from UtilFunctions import GeneFilter, FeatureNormalization



"""
============================
Method Name: Feature Cluster
============================

Method Description: This method clusters features of the input X using agglomerative hierarchical 
                    clustering using Ward linkage and returns the labels of the clustered features 
                    and the variance of the feature cluster centers.

Arguments:
========== 
X           -   (n x d) input matrix, where n is number of samples, d is the dimension
n_featclust -   the number of group to cluster the features into


Returns:
========
feat_labels -   d-dimensional vector containing the cluster labels of the features. 
                Labels are from 0 to n_featclust-1
feat_var    -   vector of size n_featclust containing the variance of the cluster centers 

"""

def FeatureCluster(X, n_featclust):

    n_samples, _ = X.shape
    fa = FeatureAgglomeration(n_clusters = n_featclust, affinity = "euclidean", linkage="ward")
    fa.fit(X)
    feat_labels = fa.labels_
    feat_means = np.zeros([n_samples, n_featclust])
    feat_var = np.zeros(n_featclust)
        
    for i in range(n_featclust):
        mask = feat_labels == i
        chunk = X[:, mask]
        feat_means[:, i] = np.mean(chunk, axis = 1) # Cluster center computed by computing the mean. 
        feat_var[i] = np.var(feat_means[:, i])

    return feat_labels, feat_var


"""
============================
Method Name: FeatClust
============================

Method Description: This method implements the proposed FeatClust approach. 

Arguments:
========== 
sc_obj              -   A single cell object which contains the data and metadata of genes and cells
n_clusters          -   The number of group to cluster the data into. This is in the range 2 <= n_clusters < n. Default (10)
gene_filt           -   A boolean variable. Gene filtering is done if True. Default (True)
apply_normalization -   A boolean variable. Data is normalized if True. Default (True)
normalization       -   A string variable. The type of normalization to apply. Options are "l2", "mean" and "norm6". Default ("l2")
k                   -   The number of groups to divide the features of the dataset into. Valid range n_clusters <= k < d. DefaultDefault (10)


Returns:
========
sc_obj              -   The single cell object containing the cluster labels in the CellData assay. A column is added to the cell data 
                        assay containing the cluster labels. The column name is the method name 'FeatClust'. 

"""

def FeatClust(  sc_obj,                                
                n_clusters = 10, 
                gene_filt = True, 
                apply_normalization = True, 
                normalization = "l2", 
                k = 10):

    method_name = "FeatClust"

    # Check and remove columns added to genedata and celldata by this method:
    sc_obj.removeCellDataColumn(method_name)
    X = sc_obj.getLogCounts()

    # Gene filtering
    if (gene_filt == True):    
        X_filt = GeneFilter(X, pct_dropout_min = 10, pct_dropout_max = 90)
    else:
        X_filt = X
    
    # Normalization
    if (apply_normalization == True):
        X_nrm = FeatureNormalization(X_filt, normalization)
    else:
        X_nrm = X_filt
    
    X_nrm = X_nrm.T # (samples, features)

    # Number of samples and features after gene filtering and normalization
    n_samples, n_features = X_nrm.shape
    
    n = np.arange(n_clusters, k+1)
    n = np.flip(n)

    s_score = np.zeros(n.shape[0])
    pred_labels = np.zeros([n_samples, n.shape[0]])

    model = AgglomerativeClustering(n_clusters=n_clusters)

    print("Computing clusters . . .")
    i = 0

    feat_labels, feat_var = FeatureCluster(X_nrm, k)
    idx = np.argsort(feat_var)
    mask = np.ones(n_features, dtype = bool)

    for j in n:

        mask = (feat_labels != idx[i]) & mask
        X_red = X_nrm[:, mask]

        pred_labels[:, i] = model.fit_predict(X_red)
        s_score[i] = silhouette_score(X_filt.T, pred_labels[:, i])

        i = i + 1

    # Select the label with the highest average silhouette coefficient
    index = np.argmax(s_score)

    # Save the final result
    print("Saving final cluster labels in Single Cell object . . .")
    sc_obj.addCellDataColumn(col_data = pred_labels[:, index], col_name = method_name)
    
    return sc_obj

