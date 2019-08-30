# Import core python libraries
import pandas as pd
import numpy as np

# Import our methods
from featclust import Cluster, PCAPlot, SilhouettePlot, MultiSilhouettePlot, GeneFilter, tSNEPlot
from singlecell import SingleCell

# Import function for metric computation
from sklearn.metrics.cluster import adjusted_rand_score

# Dataset parameters as python dictionary
dataset_params = {  'biase': {'label': 'cell_type1', 'nc': 3}}

dataset = 'biase'

data_path = "data/" + dataset + '/' + dataset + "_data.csv"
celldata_path = "data/" + dataset + '/' + dataset + "_celldata.csv"
genedata_path = "data/" + dataset + '/' + dataset + "_genedata.csv"

data = pd.read_csv(data_path, index_col=0)
celldata = pd.read_csv(celldata_path, index_col=0)
genedata = pd.read_csv(genedata_path, index_col = 0)
        
# Create a single cell object
sc = SingleCell(dataset,data,celldata,genedata)

# Get parameters from the dictionary
label = dataset_params[dataset]['label']
nc = dataset_params[dataset]['nc']

# Select only the first 49 cells from the biase dataset
if (dataset == "biase"):
    sc = sc[0:49]


# Perform Gene Filtering
sc = GeneFilter(sc, min_cells = int((10/100) * sc.dim[1]), max_cells = int((90/100) * sc.dim[1]))

# Perform clustering 
Cluster  (  sc,                                
            n_clusters = [2, 3, 4, 5], 
            apply_normalization = True, 
            normalization = "l2", 
            k = int((20/100) * sc.dim[1]))

# Compute and print the Adjusted Rand Score
print ("Adjusted Rand Index: ")
print(adjusted_rand_score(sc.getNumericCellLabels(label), sc.getCellColumnValues("FeatClust_" + str(nc) + "_Clusters")))


# Plot dataset with original cluster labels
#plot(tSNEPlot(sc, color_by = label, dist_metric="linear"), filename="./graphics/original-pcaplot.html")

# Do a Scatter Plot of the results of clustering
#plot(PCAPlot(sc, color_by = "FeatClust_" + str(nc) + "_Clusters", dist_metric="linear"), filename="./graphics/FeatClust-pcaplot.html")

# Plot silhouette and scatter plots as subplots and save plot in pdf
#fig = SilhouettePlot(sc, cluster_label="FeatClust_" + str(nc) + "_Clusters")
#plot(fig, filename="./graphics/silhouette_plot.html")
#write_image(fig, './graphics/silhouette-scatter_plot.pdf')
labs = ["FeatClust_" + str(2) + "_Clusters", "FeatClust_" + str(3) + "_Clusters", "FeatClust_" + str(4) + "_Clusters", "FeatClust_" + str(5) + "_Clusters"]
fig = MultiSilhouettePlot(sc, cluster_labels = labs)
fig.show()
