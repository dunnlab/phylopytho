# phylopytho

Phylogenetic tools written in python by the [Dunn Lab](http://dunnlab.org/).


To prepare a conda environment (before first use):

    conda create -n phylopytho -c bioconda dendropy pytest

To activate the conda environment (before each use):

    conda activate phylopytho


## treeprune

treeprune is a method and tool to decompose gene trees into subtrees that have no more than one tip per species. It does this by traversing the tree, removing the maximally inclusive subtree that has no more than one tip per species, and then repeating this process on the remaining tree until it too has no more than one tip per species. If there are two subtrees with no more than one tip per species and they both have the same size and it is the maximum, they are both pruned.

treeprune is not sensitive to the position of the tree root.

Because treeprune does not require any information beyond the gene trees, such as a species tree, it is appropriate for processing gene trees in preparation for construction of supermatrices for phylogenomic inference.

Gene trees often have clades of multiple sequences from the same species. These are often technical artifacts, especially in trees based on transcriptome assemblies that can include splice variants per gene. These monophyletic groups from the same species reduce the size of the maximally inclusive subtrees, and lead treeprune to prune an input tree into a large number of small trees. treeprune can therefore be preceded by monophyly masking. This method identifies all clades that contain multiple tips exclusively from one species, and randomly masks all but one of those tips.

Please cite the following two papers when publishing analyses that use treeprune. The first manuscript is the original description of the method, the second is the implementation of the method provided in our tool [Agalma](https://bitbucket.org/caseywdunn/agalma/src/master/), from which the stand-alone code provided here is derived.

- Hejnol, A, M Obst, A Stamatakis, M Ott, G Rouse, G Edgecombe, P Martinez, J Baguñà, X Bailly, U Jondelius, M Wiens, WEG Müller, Elaine Seaver, WC Wheeler, MQ Martindale, G Giribet, and CW Dunn (2009) Assessing the root of bilaterian animals with scalable phylogenomic methods. Proc. R. Soc. B. 276:4261-4270 http://dx.doi.org/10.1098/rspb.2009.0896

- CW Dunn, M Howison, and F Zapata (2013) Agalma: an automated phylogenomics workflow. BMC Bioinformatics 14:330. http://dx.doi.org/10.1186/1471-2105-14-330

Here is the original description from that first paper:

> In the first step, monophyly masking, all but one sequence were deleted in clades of sequences derived from the same taxon. The retained sequence was chosen at random. Paralogue pruning, the next step, consisted of identifying the maximally inclusive subtree that has no more than one sequence per taxon. This tree is then pruned away for further analysis, and the remaining tree is used as a substrate for another round of pruning. The process is repeated until the remaining tree has no more than one sequence per taxon. If there were multiple maximally inclusive subtrees of the same size in a given round, then they were all pruned away at the same time.


The treeprune code presented here was derived directly from the Agalma code, and includes the following modifications:
- Code was update from python 2 to python 3
- Only the code and dependencies relevant to treeprune were retained, greatly reducing dependencies
- Code was added so that treeprune can be used as a standalone tool rather than as part of the agalma pipeline.

### Example of command line treeprune analysis

In the `src/phylopytho` directory within this repository, first activate the conda environment as described above, and then run:

    python treeprune.py data/gene_trees.tre pruned_trees.tre

This will generate a set of pruned trees from an example set of trees included in the module.

To get a full description of options run:

    python treeprune.py -h

## Development

### Running tests

In the root of this repository, run:

    python -m pytest
