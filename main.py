# Import core python libraries
import pandas as pd
import numpy as np

# Import our methods
from FeatClust import FeatClust
from SingleCell import SingleCell

# Import function for metric computation
from sklearn.metrics.cluster import adjusted_rand_score

# Import Libraries for Plotting Functions
from PlotFunctions import PCAPlot, SilhouetteAndScatterPlot
from plotly.offline import plot
from plotly.io import write_image

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

# Exclude the following cells from analysis of biase dataset
if (dataset == "biase"):
    cells = ["blast"]
    sc.setExcludedCells(cells, label)


# Perform clustering 
FeatClust(  sc,                                
            n_clusters = nc, 
            gene_filt = True, 
            apply_normalization = True, 
            normalization = "l2", 
            k = int((20/100) * sc.dim[1]))

# Compute and print the Adjusted Rand Score
print ("Adjusted Rand Index: ")
print(adjusted_rand_score(sc.getNumericCellLabels(label), sc.getCellColumnValues("FeatClust")))


# Plot dataset with original cluster labels
plot(PCAPlot(sc, color_by = label, dist_or_kernel="linear"), filename="original-pcaplot.html")

# Do a Scatter Plot of the results of clustering
plot(PCAPlot(sc, color_by = "FeatClust", dist_or_kernel="linear"), filename="FeatClust-pcaplot.html")

# Plot silhouette and scatter plots as subplots and save plot in pdf
fig = SilhouetteAndScatterPlot(sc, n_clusters=nc, method_name="FeatClust")
plot(fig, filename="silhouette-scatter_plot.html")
write_image(fig, 'silhouette-scatter_plot.pdf')
