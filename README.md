# FeatClust
FeatClust: Clustering of Small Sample Single-Cell RNA-Seq datasets via Feature Clustering and Selection

### Requirements
You need to have latest update of python 3 and the following python libraries installed:
1. NumPy 1.15.4
2. Pandas 0.24.0
3. scikit-learn 0.20.2
4. scipy 1.2.0
5. plotly 3.6.0

### Installation
Simply clone or download the repository to your computer. 

### Usage
Navigate to the folder where you downloaded FeatClust and run the main.py python script.

### File Description
The FeatClust package contains several python files. A description of the files if as follows:
1. SingleCell.py - This file contains the implementation of SingleCell class which manages single cell datasets. 
2. PlotFunctions.py - This file contains functions for creating two types of plots. 
3. UtilFunctions.py - Contains some utility functions for pre-processing and normalization of single cell dataset.
4. PairwiseDists.py - Contains functions for computing pairwise distances and kernels between cells. Various distances are included.
5. Transformations.py - Contains functions to perform PCA (and in the future other types of transformations) on the dataset. 
6. FeatClust.py - Contains the implementation of the FeatClust method. 
