# Import Packages
import pandas as pd
import numpy as np 

from sklearn.cluster import FeatureAgglomeration, AgglomerativeClustering


from sklearn.metrics import silhouette_samples, silhouette_score

from UtilFunctions import GeneFilter, FeatureNormalization



def FeatureScores(X, labels, n_clusters):

    s_score = silhouette_samples(X, labels)
    s_average = np.mean(s_score)

    n_groups = 0                                # Number of groups in which the silhouette score for any sample is > average for all samples 
    feat_score = np.zeros(n_clusters)

    for i in range(n_clusters):
        mask1 = labels == i
        s_score_i = s_score[mask1]
        mask2 = s_score_i > s_average

        feat_score[i] = np.sum(mask2)/np.sum(mask1)
        # if (np.amax(s_score_i) > s_average):    # If the silhouette score of any sample is > average for all samples
        #     n_groups = n_groups + 1
        
    # return (n_groups/n_clusters)       # return both scores
    if (s_average < 0):
        return 0
    else:
        return s_average



def FeatureExtraction(X, n_featclust):

    """
    X           - (n x d) input matrix, where n is number of samples, d is the dimension
    n_featclust - number of groupings of features in data

    *check if any features feature groups dont have any features 
    """
    n_samples, _ = X.shape
    fa = FeatureAgglomeration(n_clusters = n_featclust, affinity = "euclidean", linkage="ward", memory="featagg-tree", compute_full_tree = True)
    #fa = FeatureAgglomeration(n_clusters = n_clusters, affinity = "euclidean", linkage="ward")
    fa.fit(X)
    feat_labels = fa.labels_
    feat_means = np.zeros([n_samples, n_featclust])
    #X_red = np.zeros([n_samples, n_featcust])
    feat_var = np.zeros(n_featclust)
        
    for i in range(n_featclust):
        mask = feat_labels == i
        chunk = X[:, mask]
        feat_means[:, i] = np.mean(chunk, axis = 1)
        feat_var[i] = np.var(np.mean(chunk, axis = 1))

    idx = np.argmin(feat_var)
    # print("variance =", feat_var[idx])
    mask2 = feat_labels != idx
    X_red = X[:, mask2]

    return X_red
    #return feat_means, feat_var, feat_labels

# Feature Agglomeration Hierarcihical Clustering
def FeatClust(  sc_obj,                                
                n_clusters = 10, 
                gene_filt = True, 
                apply_normalization = True, 
                normalization = "l2", 
                compare_with = 'cell_type1',
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
    
    # print(X_nrm[0,0])
    X_nrm = X_nrm.T # (samples, features)

    # Number of samples and features after gene filtering and normalization
    n_samples, _ = X_nrm.shape
    
    # r = int((top_n_pct/100) * n_samples)
    # print("No of features through %: ", r)
    n = np.arange(n_clusters, k+1)
    n = np.flip(n)

    s_score = np.zeros(n.shape[0])
    pred_labels = np.zeros([n_samples, n.shape[0]])


    # means, var, labs = FeatureExtraction(X_nrm, k)

    model = AgglomerativeClustering(n_clusters=n_clusters)


    print("Computing clusters . . .")
    k = 0

    for j in n:

        X_red = FeatureExtraction(X_nrm, j)

        pred_labels[:, k] = model.fit_predict(X_red)
        s_score[k] = FeatureScores(X_filt.T, pred_labels[:, k], n_clusters)
        X_nrm = X_red
        k = k + 1

    # Select the label with the highest average silhouette coefficient
    index = np.argmax(s_score)

    # Save the final result
    print("Saving final cluster labels in Single Cell object . . .")
    sc_obj.addCellDataColumn(col_data = pred_labels[:, index], col_name = method_name)
    
    return sc_obj

