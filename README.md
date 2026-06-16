# VitisViniferaPhytree
VitisViniferaPhytree computes the phylogenetic distance and builds the phylogenetic trees of the dataset [DATASET](https://doi.org/10.1007/978-3-031-64636-2_10).

# Requirements
The following tools are required:
* MATLAB with the following toolboxes:
  * Bioinformatic Toolbox
  * Statistics and Machine Learning Toolbox
* R with the following library:
  * ape

# Functions Reference

Each step of the pipeline is isolated in a MATLAB function file. Below is a brief overview of each file:
* **distAlleles.m** computes the distance between the given alleles.
* **distSingleGene.m** computes the distance between taxa of a specific gene.
* **distance.m** computes a loop for each chromosome, where it computes the distances of all the chromosome's genes and the distance for the single chromosome. It also computes the distance between taxa considering all chromosomes.
* **meanMatrix.m** computes the arithmetic mean of the matrix that are included in a folder.
* **phytreeBuildAndAnalysis.m** builds the phylogenetic trees and computes topology and evolution distance on these trees.
* **createFolder.m** creates the folder where the distances are saved.

# Trees
The trees in the folders treeNJ are built using Neighbor Joining methods, and the ones in treeME using BME.

# Matrix
Matrix folder contains the analysis performed on the trees:
* **NJTopology.txt** : Contains the topology distance (Generalized Robinson-Foulds distance) between all trees built with NJ method.
* **NJEvolution.txt**: For each phylogenetic tree, an evolutionary distance matrix is computed. Then for each possibile pair of trees the difference between their related evolutionary matrices is calculated. The Frobenius norm of each difference matrix is then computed and these values are the entries for this matrix.
* **NJvsBMETopology**: Contains the topology distances between corresponding trees built with the two different methods (NJ vs. BME).
* **NJvsBMEEvolution.txt**: Follows the same logic as 'NJEvolution.txt',  but the analysed pairs are corresponding trees built with different reconstruction methods.
