from SingleCell import SingleCell
from Transformations import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score, silhouette_samples

import numpy as np
import plotly.graph_objs as go
from plotly import tools

import colorlover as cl

from UtilFunctions import GeneFilter, FeatureNormalization

def PCAPlot(sc_obj, color_by, dist_or_kernel = "linear"):

    X = sc_obj.getLogCounts() # X is (features, samples)
    X = GeneFilter(X, pct_dropout_min = 10, pct_dropout_max = 90)
    X = FeatureNormalization(X, norm = "mean")
    X_red = PCA(X = X, dist_or_kernel = dist_or_kernel, n_comp = 2)
    x = X_red[0, :]
    y = X_red[1, :]

    colors = cl.scales['8']['qual']['Set1']

    # To Do: Check if color_by is valid
    data = []
    cell_types = sc_obj.getDistinctCellTypes(color_by)
    cell_labels = sc_obj.getNumericCellLabels(color_by)

    for i in range(cell_types.size):
        mask = (cell_labels == i)
        trace = go.Scatter(x = x[mask],
                           y = y[mask], 
                           mode = 'markers', 
                           showlegend = True, 
                           name = cell_types[i],
                           marker = dict(size = 10, color = colors[i]))
        data.append(trace)
        
    layout = go.Layout(title = "PCA Plot for Dataset: " + sc_obj.dataset, 
                       xaxis = dict(title = 'PC1',
                                    zeroline = False, 
                                    showgrid = True),
                       yaxis = dict(title = 'PC2',  
                                    zeroline = False, 
                                    showgrid = True))
    
    fig = go.Figure(data = data, layout = layout)
    
    return fig
    
"""
def tSNEPlot(sc_obj, color_by):

    sc_obj.setLabelColumn(color_by)
    tsne = TSNE(n_components = 2)
    X_red_tsne = tsne.fit_transform(sc_obj.getLogCounts().T)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for i in range(sc_obj.cell_types.size):
        mask = sc_obj.celldata["num_cell_labels"].values == i
        ax.scatter(X_red_tsne[mask, 0], X_red_tsne[mask, 1], label = sc_obj.cell_types[i])
    
    ax.set_title("tSNE Plot for Dataset: " + sc_obj.dataset)
    ax.set_xlabel("Dim1")
    ax.set_ylabel("Dim2")
    ax.legend(loc = 0)
    
    xevents = EventCollection(X_red_tsne[:, 0], linelength=0.01)
    yevents = EventCollection(X_red_tsne[:, 1], linelength=0.01, orientation="vertical")
    ax.add_collection(xevents)
    ax.add_collection(yevents)
    
    plt.show()


     
                                    tickmode = 'array', 
                                    ticks = 'inside', 
                                    ticklen = 10,
                                    tickvals = x, 
                                    showticklabels = False, 
"""


def SilhouetteAndScatterPlot(sc_obj, n_clusters, method_name, pca_dist_or_kernel = "linear"):

    X = sc_obj.getLogCounts() # X is (features, samples)
    X_filt = GeneFilter(X, pct_dropout_min = 10, pct_dropout_max = 90)
    X_nrm = FeatureNormalization(X_filt, norm = "mean")
    X_red = PCA(X = X_nrm, dist_or_kernel = pca_dist_or_kernel, n_comp = 2)
    x = X_red[0, :]
    y = X_red[1, :]

    colors = cl.scales['8']['qual']['Set1']
    labels = sc_obj.getCellColumnValues(method_name)

    # Create a subplot with 1 row and 2 columns
    fig = tools.make_subplots(rows=1, cols=2)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    fig['layout']['xaxis1'].update(title='silhouette coefficient', zeroline = True)
    fig['layout']['yaxis1'].update(title='samples',showticklabels=True, zeroline = True, range=[0, len(X_nrm.T)])
    fig['layout']['xaxis2'].update(title='PC1',zeroline=False)
    fig['layout']['yaxis2'].update(title='PC2',zeroline=False)
    # fig['layout'].update(title="Silhouette Analysis for Method '" + method_name + "' on Dataset '" + sc_obj.dataset + "'")

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    silhouette_avg = silhouette_score(X_filt.T, labels)
    print("For n_clusters =", n_clusters, "The average silhouette_score is :", silhouette_avg)

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X_filt.T, labels)
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
                                 line=dict(width=0.5, color = colors[i]),
                                 fill='tozerox')

        fig.append_trace(filled_area, 1, 1)
        
        # Compute the new y_lower for next plot
        y_lower = y_upper + 0  # 10 for the 0 samples
        

    # The vertical line for average silhouette score of all the values
    axis_line = go.Scatter(x=[silhouette_avg, silhouette_avg],
                           y=[0, len(X_nrm.T)],
                           showlegend=False,
                           mode='lines',
                           line=dict(color="red", dash='dash',width =2) )

    fig.append_trace(axis_line, 1, 1)
    
    # 2nd Plot showing the actual clusters formed

    for i in range(n_clusters):
        mask = (labels == i)
        trace = go.Scatter(x = x[mask],
                           y = y[mask], 
                           mode = 'markers', 
                           showlegend = True, 
                           name = "cell_type-" + str(i+1),
                           marker = dict(size = 5, color = colors[i]))
        fig.append_trace(trace, 1, 2)

    return fig