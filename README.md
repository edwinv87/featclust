# FeatClust

FeatClust: Clustering of Small Sample Single-Cell RNA-Seq datasets via Feature Clustering and Selection.
FeatClust is a python package for clustering and visualization of single-cell RNA-Seq datasets. It is implemented using popular python libraries such as numpy, pandas and scikit-learn. FeatCust can also be used to visualize high dimensional data in 2D using PCA and tSNE and generate silhouette plots to determine the appropriate number of clusters in a dataset visually.  

## Requirements

FeatClust uses the following python packages:

1. numpy
2. pandas
3. scikit-learn
4. scipy
5. plotly
6. singlecell

## Installation

FeatClust is available on Anaconda cloud <https://anaconda.org/edwinvans/featclust> and can be easily installed using the Anaconda prompt if you are using Anaconda distribution. On the Anaconda prompt change to the environment in which you want to install FeatClust by typing `conda activate env`, where env is the environment. Then  type the following command to install FeatCust and all its dependencies:

`conda install -c edwinvans featclust`

Alternatively, you can also add the edwinvans channel in Anaconda distribution using the following code so that package updates can be downloaded and installed automatically.

`conda config --add channels edwinvans`

## Usage

A Python 3 Jupyter Notebook explaining hot to use FeatClust is available at: <https://anaconda.org/edwinvans/installing_and_using_featclust/notebook>

## Corresponding Paper

The corresponding paper of the FeatClust method is available at <https://link.springer.com/chapter/10.1007/978-3-030-29894-4_36>

## Contact

For reporting bugs or help regarding any aspect of the FeatClust method, please email: vans.edw@gmail.com
