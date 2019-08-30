# Import libraries
import pandas as pd
import numpy as np
from singlecell import SingleCell

from sklearn.metrics.cluster import adjusted_rand_score

from UtilFunctions import GeneFilter, FeatureNormalization
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.feature_selection import SelectKBest, f_classif, mutual_info_classif


dataset_params = {  'biase': {'label': 'cell_type1', 'nc': 3},
                    'yan': {'label': 'cell_type2', 'nc': 7},
                    'goolam': {'label': 'cell_type1', 'nc': 5},
                    'fan': {'label': 'cell_type1', 'nc': 6},
                    'treutlein': {'label': 'cell_type1', 'nc': 5},
                    'ting': {'label': 'cell_type2', 'nc': 7}}


# Load Data
dataset = 'yan'
data_path = "../../New_Datasets/" + dataset + '/' + dataset + "_data.csv"
celldata_path = "../../New_Datasets/" + dataset + '/' + dataset + "_celldata.csv"
genedata_path = "../../New_Datasets/" + dataset + '/' + dataset + "_genedata.csv"

data = pd.read_csv(data_path, index_col=0)
celldata = pd.read_csv(celldata_path, index_col=0)
genedata = pd.read_csv(genedata_path, index_col = 0)
    
# Create a single cell object
sc = SingleCell(dataset,data,celldata,genedata)

# Get parameters from the dictionary
label = dataset_params[dataset]['label']
nc = dataset_params[dataset]['nc']

# Exclude the following cells from analysis of biase dataset
if (dataset == "biase"):
    cells = ["blast"]
    sc.setExcludedCells(cells, label)

print("Dimension of data: ", sc.dim)



# Parameters
gene_filt = True
apply_normalization = True,
normalization = "l2" 
classification_type = "f_classif"
top_n_pct = 10


method_name = "AnovaKMeans"

# Check and remove columns added to genedata and celldata by this method:
sc.removeCellDataColumn(method_name)
X = sc.getLogCounts()

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

# FIND TEMPORARY CLASS LABELS
# Using Hierarchical Clustering
print("Computing Temporary Clusters . . .")
model_hc = AgglomerativeClustering(n_clusters = nc)
temp_labels = model_hc.fit_predict(X_nrm)
temp_s_score = silhouette_score(X_filt.T, temp_labels)
    
print("Temp. Label ARI: ", adjusted_rand_score(sc.getNumericCellLabels(label), temp_labels))
print("Temp. Label S-Score: ", temp_s_score)

print("Performing Feature Selection . . .")
if (classification_type == "f_classif"):
    feat_sel = SelectKBest(f_classif, k="all")
    
elif (classification_type == "mutual_info_classif"):
    feat_sel = SelectKBest(mutual_info_classif, k="all")

feat_sel.fit(X_nrm, temp_labels)
feature_scores = feat_sel.scores_
idx = np.argsort(feature_scores, kind='mergesort')
idx = idx[::-1]

r = int((top_n_pct/100) * X_nrm.shape[0])
# print("No of features through %: ", r)
n = np.arange(1, r+1)
med_ari = np.zeros(n.shape[0])
s_score = np.zeros(n.shape[0])
pred_labels = np.zeros([X_nrm.shape[0], n.shape[0]])

model_km = KMeans(n_clusters=nc, n_init=100, max_iter=10000, n_jobs=-1)

print("Computing clusters using best n features . . .")
k = 0

for j in n:
    X_red = X_nrm[:, idx[0:j]]
    pred_labels[:, k] = model_km.fit_predict(X_red)
    s_score[k] = silhouette_score(X_filt.T, pred_labels[:, k])
    med_ari[k] = adjusted_rand_score(sc.getNumericCellLabels(label), pred_labels[:, k])
    print("N = ", j, "S-SCORE = ", s_score[k], "ARI = ", med_ari[k])
    k = k + 1

# Select the label with the highest average silhouette coefficient
index = np.argmax(s_score)

# Save the final result
print("Saving final cluster labels in Single Cell object . . .")
sc.addCellDataColumn(col_data = pred_labels[:, index], col_name = method_name)

import plotly.plotly as py
import plotly.graph_objs as go
from plotly import tools

# Silhouette Analysis
def SilhouettePlot(X, labels, n_clusters):

     # Create a subplot with 1 row and 2 columns
    fig = tools.make_subplots(rows=1, cols=2,
                              subplot_titles=('Silhouette Plot for Various Clusters',
                                              'Visualization of Clustered Data'))

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    fig['layout']['xaxis1'].update(title='silhouette coefficient',range=[-0.1, 1], zeroline = True)
    fig['layout']['yaxis1'].update(title='Samples',showticklabels=True,zeroline = True, range=[0, len(X)])
    fig['layout']['xaxis2'].update(title='PC1',zeroline=False)
    fig['layout']['yaxis2'].update(title='PC2',zeroline=False, range=[0, 5])
    fig['layout'].update(title="Silhouette analysis for KMeans clustering on sample data with n_clusters = %d" % n_clusters)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    silhouette_avg = silhouette_score(X, labels)
    print("For n_clusters =", n_clusters, "The average silhouette_score is :", silhouette_avg)

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, labels)
    y_lower = 0

    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[labels == i]
        ith_cluster_silhouette_values.sort()
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
        
        filled_area = go.Scatter(y=np.arange(y_lower, y_upper),
                                 x=ith_cluster_silhouette_values,
                                 mode='lines',
                                 showlegend=False,
                                 line=dict(width=0.5),
                                 fill='tozerox')

        fig.append_trace(filled_area, 1, 1)
        
        # Compute the new y_lower for next plot
        y_lower = y_upper + 0  # 10 for the 0 samples
        

    # The vertical line for average silhouette score of all the values
    axis_line = go.Scatter(x=[silhouette_avg, silhouette_avg],
                           y=[0, len(X)],
                           showlegend=False,
                           mode='lines',
                           line=dict(color="red", dash='dash',width =2) )

    fig.append_trace(axis_line, 1, 1)
    
    # 2nd Plot showing the actual clusters formed
    clusters = go.Scatter(x=X[:, 0], 
                          y=X[:, 1], 
                          showlegend=False,
                          mode='markers',
                          marker=dict(size=8))
    fig.append_trace(clusters, 1, 2)
    
    """
    # Labeling the clusters
    centers_ = clusterer.cluster_centers_
    # Draw white circles at cluster centers
    centers = go.Scatter(x=centers_[:, 0], 
                         y=centers_[:, 1],
                         showlegend=False,
                         mode='markers',
                         marker=dict(color='green', size=10,
                                     line=dict(color='black',
                                                             width=1))
                        )

    fig.append_trace(centers, 1, 2)
    """



    return fig

from plotly.offline import plot

plot(SilhouettePlot(X_filt.T, pred_labels[:, 3], nc))