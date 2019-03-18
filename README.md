# FeatClust
FeatClust: Clustering of Small Sample Single-Cell RNA-Seq datasets via Feature Clustering and Selection.
FeatClust is for clustering and visualization of SingleCell RNA-Seq datasets. It is implemented in python programming language. 

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
Navigate to the folder where you downloaded FeatClust and run the main.py python script. Note a small dataset has been provided here for testing. To use your own datset, you may need to change the path to point to your dataset. Note that FeatClust will open the plots in a browser environment and the same plots will be saved in your folder as 'html' file to be accessed later.  

### File Description
The FeatClust package contains several python files. A description of the files is as follows:
1. SingleCell.py - This file contains the implementation of SingleCell class which manages single cell datasets. 
2. PlotFunctions.py - This file contains functions for creating two types of plots. 
3. UtilFunctions.py - Contains some utility functions for pre-processing and normalization of single cell dataset.
4. PairwiseDists.py - Contains functions for computing pairwise distances and kernels between cells. Various distances/kernels are included.
5. Transformations.py - Contains functions to perform PCA (and in the future other types of transformations) on the dataset. 
6. FeatClust.py - Contains the implementation of the FeatClust method. 

### Contact
For reporting bugs or help regarding any aspect of the FeatClust method, please email: vans.edw@gmail.com 
