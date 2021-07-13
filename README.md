# phylopytho

Phylogenetic tools written in python.


To prepare a conda environment:

    conda create -n phylopytho -c bioconda biopython dendropy pytest
    conda activate phylopytho


## treeprune

treeprune is a method and tool to decompose gene trees into subtrees that have no more than one tip per species. It does this by traversing the tree, removing the maximally inclusive subtree that has no more than one tip per species, and then repeating this process on the remaining tree until it too has no more than one tip per species. If there are two subtrees with no more than one ip per species and they both have the same size and it is the maximum, they are both pruned.


treeprune is not sensitive to the position of the tree root.

Because treeprune does not require any information beyone the gene trees, such as a species tree, it is appropriate for processing gene trees in preparation for construction of supermatrices for phylogenomic inference.

Please cite the following papers when publishing analyses that use treeprune:




## Running tests

In the root of this repository, run:

    python -m pytest


