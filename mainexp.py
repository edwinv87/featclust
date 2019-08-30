# Import libraries
import pandas as pd
import numpy as np

from FeatClust import FeatClust
from singlecell import SingleCell

from sklearn.metrics.cluster import adjusted_rand_score

dataset_params = {  'biase': {'label': 'cell_type1', 'nc': 3},
                    'yan': {'label': 'cell_type2', 'nc': 7},
                    'goolam': {'label': 'cell_type1', 'nc': 5},
                    'fan': {'label': 'cell_type1', 'nc': 6},
                    'treutlein': {'label': 'cell_type1', 'nc': 5},
                    'ting': {'label': 'cell_type2', 'nc': 7}}

n_trials = 1
datasets = ['biase', 'yan','treutlein', 'goolam', 'fan', 'ting']
ari = np.zeros([n_trials, len(datasets)])

col = 0
for dataset in datasets:
    # Load Data
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

# Select only the first 49 cells from the biase dataset
if (dataset == "biase"):
    sc = sc[0:49]

    print("Dimension of data: ", sc.dim)

    row = 0

    for i in range(n_trials):
    # Perform clustering 
        FeatClust(  sc,                                
                    n_clusters = nc, 
                    gene_filt = True, 
                    apply_normalization = True, 
                    normalization = "l2", 
                    k = int((80/100) * sc.dim[1]))

        # Compute and print the Adjusted Rand Score
        ari[row, col] = adjusted_rand_score(sc.getNumericCellLabels(label), sc.getCellColumnValues("FeatClust"))
        print ("Adjusted Rand Index: ")
        print(ari[row, col])
        row = row + 1

    col = col + 1


df = pd.DataFrame(data = ari, columns = datasets)
df.to_excel("results/results-featclust_80.xlsx")

