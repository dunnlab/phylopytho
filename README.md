# phylopytho

Phylogenetic tools written in python.


To prepare a conda environment:

    conda create -n phylopytho -c bioconda dendropy pytest
    conda activate  phylopytho


## treeprune

treeprune is a method and tool to decompose gene trees into subtrees that have no more than one tip per species. It does this by traverseeing the tree, removing the maximally inclusive subtree that has no more than one tip per species, and then repeating this process on the remaining tree until it too has no more than one taxon per species.

treeprune is not sensitive to the position of the tree root.

Because it does not require any information beyone the gene trees, such as a species tree, it is appropriate for processing gene trees in preparation for construction of supermatrices for phylogenomic inference.


## Running tests

In the root of this repository, run:

    pytest


