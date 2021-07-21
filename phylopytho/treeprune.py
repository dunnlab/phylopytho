#!/usr/bin/env python

import dendropy
import argparse
from collections import OrderedDict

# import sample data, could be used for an --example option
# import pkg_resources
# example_tree_file = pkg_resources.resource_string(__name__, "data/gene_trees.tre")

def get_species(taxon):
	"""
	Parses the species from a taxon name in the format 'species@sequence_id'.
	"""
	return str(taxon).partition('@')[0]


def get_taxa(nodes):
	"""
	Returns a list of taxa names given a list of nodes.
	"""
	return [node.taxon.label for node in nodes]

def validate_names(tree):
	"""
	Returns a list of any tip names that do not conform to the
	taxon@id format. An empty list indicates that all names conform.
	"""
	invalid_names = list()
	for tip_name in get_taxa(tree.leaf_nodes()):
		valid = True
		fields = tip_name.split("@")
		if len(fields) != 2:
			valid = False
		else:
			if not len(fields[0]) > 0:
				valid = False
			if not len(fields[0]) > 0:
				valid = False
		if not valid:
			invalid_names.append(tip_name)

	return invalid_names

def count_species(taxa):
	"""
	Returns the count of species given a list of taxon names.
	"""
	species = set(map(get_species, taxa))
	return len(species)


def count_orthologs(taxa):
	"""
	Returns the count of orthologs given a list of taxon names. If there are as
	many species as there are taxon names, then the count is the number of
	taxon names. If there are less species than taxon names, then the count is
	0 since there must be some paralogs.
	"""
	species = list(map(get_species, taxa))
	return len(species) * int(len(species) == len(set(species)))


def monophyly_prune(tree):
	"""
	Takes a tree and identifies clades that have more than one sequence per
	taxon and prunes tip at random leaving a single representative sequence per
	taxon.
	"""

	invalid_names = validate_names(tree)
	if len(invalid_names) > 0:
		raise ValueError(f"Tip name {invalid_names[0]} does not conform to taxon@id naming requirement.")


	for node in tree.internal_nodes():
		if node.parent_node:
			tree.reroot_at_node(node)
			for leaf in tree.leaf_nodes():
				sister = leaf.sister_nodes()[0]
				if get_species(leaf.taxon) == get_species(sister.taxon):
					tree.prune_taxa([leaf.taxon], update_bipartitions=False)
	return tree


def split_tree(tree, nodes):
	"""
	Takes a rooted tree and a list of non-root internal nodes,
	and returns a set of the subtrees from clipping below
	the nodes.
	"""
	all_taxa_leafset = frozenset([nd for nd in tree.leaf_node_iter()])
	leafset_of_new_trees = set()
	for target_node in nodes:
		to_remove = []
		to_add = []
		leafset1 = frozenset([nd for nd in target_node.leaf_iter()])
		proposed_leafsets = [
				leafset1,
			    frozenset(all_taxa_leafset - leafset1),
			    ]
		for existing_leafset in leafset_of_new_trees:
			for pidx, proposed_leafset in enumerate(proposed_leafsets):
				if existing_leafset.issubset(proposed_leafsets[pidx]):
					proposed_leafsets[pidx] = proposed_leafsets[pidx] - existing_leafset
				elif proposed_leafsets[pidx].issubset(existing_leafset):
					to_remove.append(existing_leafset)
					to_add.append( existing_leafset - proposed_leafsets[pidx] )
		for leafset_taxa in to_remove:
			leafset_of_new_trees.discard(leafset_taxa)
		for leafset_taxa in to_add:
			leafset_of_new_trees.add(leafset_taxa)
		for proposed_leafset in proposed_leafsets:
			leafset_of_new_trees.add(proposed_leafset)
	results = set()
	for leafset in leafset_of_new_trees:
		if not leafset:
			continue
		induced_tree = tree.extract_tree(
				node_filter_fn=lambda nd: nd in leafset)
		results.add(induced_tree)
	return results

def paralogy_prune(tree, pruned_trees):
	"""
	Takes a tree and loops through the internal nodes (except root) and gets
	the maximum number of orthologs on either side of each node.
	"""

	invalid_names = validate_names(tree)
	if len(invalid_names) > 0:
		raise ValueError(f"Tip name {invalid_names[0]} does not conform to taxon@id naming requirement.")

	# A dictionary will keep track of node (key) and maximum number of
	# orthologs at this node (value). Internal nodes are reported by dendropy
	# in prefix order, so by using an OrderedDict, we will split the tree as
	# close to the root as possible.
	counts = OrderedDict()
	all_leaves = set(get_taxa(tree.leaf_nodes()))

	# Stop at any trees that only have two leaves.
	if len(all_leaves) < 2:
		pruned_trees.append(tree)
		#print tree.as_newick_string()
		return

	for node in tree.postorder_node_iter():
		child_leaves = set(get_taxa(node.leaf_nodes()))
		other_leaves = all_leaves - child_leaves
		counts[node] = max(
						count_orthologs(child_leaves),
						count_orthologs(other_leaves))
		node.label = str(counts[node])

	#print tree.as_ascii_plot(show_internal_node_labels=True)

	# Calculate the maximum number of orthologs across all nodes.
	all_max = max(counts.values())

	# Create a list of the nodes that have this max value, where we might split
	# the tree.
	splits = [node for node in counts if counts[node] == all_max]
	assert splits

	n_species = count_species( get_taxa(tree.leaf_nodes()) )

	if ( splits[0] == tree.seed_node ) and ( n_species > 1 ):
		assert len(splits) == 1
		# The tree is pruned and only contains orthologs.
		pruned_trees.append(tree)
	else:
		assert not tree.seed_node in splits
		# Split the tree at the first non-root split node. Because the nodes
		# are in prefix order, this should be close to the root.
		subtrees = split_tree(tree, splits)
		for subtree in subtrees:
			paralogy_prune(subtree, pruned_trees)

def main():

  # https://docs.python.org/3/howto/argparse.html
  parser = argparse.ArgumentParser(description="prune one or more gene trees into maximally inclusive subtrees with no more than one tip per species")
  parser.add_argument("input", help="input tree file, containing one or more phylogenies in newick format")
  parser.add_argument("output", help="output tree file, containing all the pruned subtrees")
  parser.add_argument("-d", "--disablemonophylyprune", help="disable monophyly_prune prior to paralogy_prune",
                    action="store_true")
  parser.add_argument("-m", "--mintreesize", type=int, default=4,
                    help="the minimum size of pruned trees to retain")
  args = parser.parse_args()
  
  input_name = args.input
  print(f'Input file: {input_name}')
  input_trees = dendropy.TreeList.get(path=input_name, schema="newick")

  output_trees = dendropy.TreeList()

  if args.disablemonophylyprune:
    print('Skipping monophyly pruning')

  n_in = 0
  n_out = 0
  n_reject = 0

  for tree in input_trees:
    n_in = n_in + 1

    if not args.disablemonophylyprune:
      tree = monophyly_prune(tree)

    pruned_trees = []
    paralogy_prune(tree, pruned_trees)

    for pruned_tree in pruned_trees:
      n_tips = len(pruned_tree.leaf_nodes())
      if n_tips >= args.mintreesize:
        n_out = n_out + 1
        output_trees.append(pruned_tree)
      else:
        n_reject = n_reject + 1	  

  output_name = args.output
  print(f'Output file: {output_name}')
  output_trees.write(path=output_name, schema="newick")

  print(f'{n_in} trees read from input file')
  print(f'{n_reject} pruned trees rejected because they had too few tips')
  print(f'{n_out} pruned trees written to the output file')

if __name__ == "__main__":
	main()