# phylopytho

Phylogenetic tools written in python by the [Dunn Lab](http://dunnlab.org/).

## Installation

To prepare a conda environment and install phylopytho (before first use):

    conda create -n phylopytho -c conda-forge -c bioconda dendropy pytest
    conda activate phylopytho
    pip install git+https://github.com/dunnlab/phylopytho.git

To activate the conda environment (before each use):

    conda activate phylopytho


## Conventions

Tip names on trees follow the convention `taxon@ID`, for example `Nematostella_vectensis@301933`. `taxon` can be any name, but is typically `Genus_species`. `ID` is a unique identifier. It can be unique to the taxon, or globally unique across all taxa. Avoid spaces or special characters in the names.


## treeprune

`treeprune`, available in the phylopytho module, is a method and tool to decompose gene trees into subtrees that have no more than one tip per species. It does this in two steps - `monophyly_prune` and `paralogy_prune`.

`paralogy_prune` traverses the tree, identifying the maximally inclusive subtree that has no more than one tip per species. It then prunes away and outputs this  maximally inclusive subtree. This process is repeated until the remaining tree has no more than one tip per species, at which time the remaining tree is also added to the output. If at any step there are two subtrees with no more than one tip per species and they both have the same size and it is the maximum, they are both pruned. `paralogy_prune` is not sensitive to the position of the tree root.

Gene trees often have clades of multiple sequences from the same species. These are often technical artifacts, especially in trees based on transcriptome assemblies that can include multiple splice variants per gene which then form monophyletic clades on the gene tree. These monophyletic groups from the same species reduce the size of the maximally inclusive subtrees, and lead `paralogy_prune` to prune an input tree into a large number of small trees. `paralogy_prune` can therefore be preceded by `monophyly_prune`. This method identifies all clades that contain multiple tips exclusively from one species, and randomly masks all but one of those tips.

Please cite the following two papers when publishing analyses that use `treeprune`. The first manuscript is the original description of the method, the second is the implementation of the method provided in our tool [Agalma](https://bitbucket.org/caseywdunn/agalma/src/master/), from which the stand-alone code provided here is derived.

- Hejnol, A, M Obst, A Stamatakis, M Ott, G Rouse, G Edgecombe, P Martinez, J Baguñà, X Bailly, U Jondelius, M Wiens, WEG Müller, Elaine Seaver, WC Wheeler, MQ Martindale, G Giribet, and CW Dunn (2009) Assessing the root of bilaterian animals with scalable phylogenomic methods. Proc. R. Soc. B. 276:4261-4270 http://dx.doi.org/10.1098/rspb.2009.0896

- CW Dunn, M Howison, and F Zapata (2013) Agalma: an automated phylogenomics workflow. BMC Bioinformatics 14:330. http://dx.doi.org/10.1186/1471-2105-14-330

Here is the original description from that first paper:

> In the first step, monophyly masking, all but one sequence were deleted in clades of sequences derived from the same taxon. The retained sequence was chosen at random. Paralogue pruning, the next step, consisted of identifying the maximally inclusive subtree that has no more than one sequence per taxon. This tree is then pruned away for further analysis, and the remaining tree is used as a substrate for another round of pruning. The process is repeated until the remaining tree has no more than one sequence per taxon. If there were multiple maximally inclusive subtrees of the same size in a given round, then they were all pruned away at the same time.

Because `treeprune` does not require any information beyond the gene trees, such as a species tree, it is appropriate for processing gene trees in preparation for construction of supermatrices for phylogenomic inference.

The `treeprune` code presented here was derived directly from the Agalma code, and includes the following modifications:
- Code was update from python 2 to python 3
- Only the code and dependencies relevant to `treeprune` were retained, greatly reducing dependencies
- Code was added so that `treeprune` can be used as a standalone tool rather than as part of the agalma pipeline.
- Monophyly masking has been renamed `monophyly_prune` to be consistent with the tern `paralogy_prune`, and to reflect the fact that it is pruning and not masking.


### Example of using treeprune within python

Both `monophyly_prune` and `paralogy_prune` are available from `phylopytho.treeprune`:

    from dendropy import Tree
    from phylopytho import treeprune

    # Build the tree
    tree_text = "(((((((A@1,B@1),C@1),A@2),A@4),(B@2,B@3)),C@2),A@3);"
    t = Tree.get_from_string( tree_text, schema='newick', rooting='force-rooted' )
    t.print_plot()

    # Perform monophyly pruning
    t = treeprune.monophyly_prune(t)
    t.print_plot()

    # Perform paralogy pruning
    pruned_trees = list()
    treeprune.paralogy_prune(t, pruned_trees)
    print(f'Number of pruned subtrees: {len(pruned_trees)}')

    # Print the trees
    for pruned_tree in pruned_trees:
      # Can only print trees with more than 1 tip
      if len(pruned_tree) > 1:
        pruned_tree.print_plot()
      else:
        print("single taxon")
    
    # Print the sets of tips in each tree
    for pruned_tree in pruned_trees:
      tips = treeprune.get_taxa(pruned_tree.leaf_nodes())
      print(tips)




### Example of command line treeprune analysis

`treeprune` can also be used at the command line as a stand-alone script when the phylopytho module is installed. By default, it applies both `monophyly_prune` and `paralogy_prune`.

To give it a try, in the root of this repository first activate the conda environment as described above, and then run:

    treeprune phylopytho/data/gene_trees.tre pruned_trees.tre

This will generate a set of pruned trees from an example set of trees included in the module.

To get a full description of options run:

    treeprune --help


## Development

### Running tests

In the root of this repository, run:

    python -m pytest
