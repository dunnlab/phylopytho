#!/usr/bin/env python
#
# BioLite - Tools for processing gene sequence data and automating workflows
# Copyright (c) 2012-2017 Brown University. All rights reserved.
#
# This file is part of BioLite.
#
# BioLite is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BioLite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioLite.  If not, see <http://www.gnu.org/licenses/>.

import dendropy
import ete3
import numpy as np
import os
from Bio import SeqIO
from collections import OrderedDict, defaultdict, namedtuple
from operator import attrgetter
from biolite import diagnostics
from biolite import utils


def tree(newick):
	"""
	Create a Dendropy tree object from a `newick` string.
	"""
	return dendropy.Tree.get_from_string(newick, schema="newick")


def ascii_plot(newick):
	"""
	Print an ascii plot from a `newick` string.
	"""
	tree = dendropy.Tree.get_from_string(newick, schema="newick")
	return tree.as_ascii_plot(show_internal_node_labels=True)


def read_tree(path):
	"""
	Read the newick file at `path` and return as a string.
	"""
	return dendropy.Tree.get_from_path(path, schema="newick")


def reroot_tree(path, outgroup):
	"""
	Read the newick file at `path` and reroot it so that the list of node
	labels in `outgroup` are in the outgroup position.
	"""
	tree = dendropy.Tree.get_from_path(path, schema="newick")
	utils.info("target tree:\n\n%s" % tree.as_ascii_plot())
	mrca = tree.mrca(taxon_labels=outgroup)
	tree.reroot_at_edge(mrca.edge)
	utils.info("rerooted at outgroup:\n\n%s" % tree.as_ascii_plot())
	return tree.as_string(schema="newick")


def get_tips(tree):
	"""
	Return a list of taxa for the leaf nodes of the `tree`.
	"""
	tree.deroot()
	nodes = tree.leaf_nodes()
	return [n.taxon.label.partition('@')[-1] for n in nodes]


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
	species = map(get_species, taxa)
	return len(species) * int(len(species) == len(set(species)))


def monophyly_masking(tree):
	"""
	Takes a tree and identifies clades that have more than one sequence per
	taxon and prunes tip at random leaving a single representative sequence per
	taxon.
	"""
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
	all_max = max(counts.itervalues())

	# Create a list of the nodes that have this max value, where we might split
	# the tree.
	splits = filter(lambda node : counts[node] == all_max, counts)
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


def identify_candidate_variants(newicks, threshold):
	"""
	Annotate each node in an ete3 `tree` with the sum of the branch lengths
	for its subtree and identify subtrees with branch length < `threshold`.
	For each species in each subtree, yield a list of the tips in the subtree
	as candidates for variants of the same gene.
	"""
	hist = {}

	for newick in newicks:
		tree = ete3.Tree(newick)
		outgroup = tree.get_midpoint_outgroup()
		if not outgroup is None:
			tree.set_outgroup(outgroup)
		for node in tree.traverse(strategy="postorder"):
			if node.is_leaf():
				node.add_feature("branchlength", 0)
				node.add_feature("under", True)
			if not node.is_leaf():
				children = node.get_children()
				branchlength = children[0].get_distance(children[1]) + children[0].branchlength + children[1].branchlength
				node.add_feature("branchlength", branchlength)
				hist[branchlength] = hist.get(branchlength, 0) + 1
				if branchlength < threshold: # adds flag for if node.bl under threshold
					node.add_feature("under", True)
				else:
					node.add_feature("under", False)

		# yield candidates for subtrees with branch length < threshold
		for node in tree.traverse(strategy="levelorder", is_leaf_fn=attrgetter("under")):
			if node.branchlength != 0 and node.under == True:
				candidates = defaultdict(set)
				for leaf in node.get_leaves():
					species, _, model_id = leaf.name.partition('@')
					model_id = int(model_id)
					candidates[species].add(model_id)
				for child in node.get_children():
					child.checked = True
				for species in candidates:
					if len(candidates[species]) > 1:
						yield candidates[species]

	diagnostics.log("histogram", hist)


def supermatrix(clusters, outdir, partition_prefix='DNA', proportion=None):
	"""
	Build a supermatrix from a list of clusters (as FASTA file names).
	"""

	Gene = namedtuple("Gene", "name seqs")

	taxa = {}
	genes = []

	for cluster in clusters:
		# Build a dictionary of taxa/sequences in this cluster
		gene = Gene(utils.basename(cluster), {})
		for record in SeqIO.parse(cluster, 'fasta'):
			gene.seqs[record.id.partition('@')[0]] = str(record.seq)
		for taxon in gene.seqs:
			taxa[taxon] = taxa.get(taxon, 0) + 1
		genes.append(gene)

	# Sort from best to worst sampled genes.
	genes = sorted(genes, key=(lambda gene: len(gene.seqs)), reverse=True)

	# Sort from best to worst sampled taxa.
	taxa = sorted(taxa, key=taxa.get, reverse=True)
	diagnostics.log("taxa", taxa)

	superseq = dict([(x, []) for x in taxa])
	occupancy = np.zeros((len(genes), len(taxa)), dtype=np.bool)
	cells = np.zeros(len(genes), dtype=np.uint32)
	sites = np.zeros(len(genes), dtype=np.uint32)

	for i, gene in enumerate(genes):

		# Find the maximum sequence length
		maxlen = max(map(len, gene.seqs.itervalues()))
		sites[i] = maxlen

		# Append sequences to the super-sequence, padding out any that
		# are shorter than the longest sequence
		for j, taxon in enumerate(taxa):
			seq = gene.seqs.get(taxon, '')
			if seq:
				superseq[taxon].append(seq)
				occupancy[i,j] = True
				cells[i] += sum(1 for c in seq if c != '-')
			if len(seq) < maxlen:
				superseq[taxon].append('-' * (maxlen - len(seq)))

	cum_occupancy = np.cumsum(np.sum(occupancy, axis=1), dtype=np.float64)
	cum_occupancy /= np.cumsum(len(taxa) * np.ones(len(clusters)))

	nsites = np.sum(sites)
	cell_occupancy = float(np.sum(cells)) / (nsites * len(taxa))

	utils.info("total genes:", len(genes))
	utils.info("total sites:", nsites)
	utils.info("total occupancy: %.2f" % cum_occupancy[-1])
	utils.info("total cell occupancy: %.2f" % cell_occupancy)

	# Use argument 'proportion' to trim supermatrix: find index in
	# cum_occupancy to trim superseq per taxon
	index = len(genes)
	if proportion:
		for i, val in enumerate(cum_occupancy):
			if val < proportion:
				index = i
				break
		if index <= 0:
			utils.die("supermatrix is empty")
		nsites = np.sum(sites[:index])
		cell_occupancy = float(np.sum(cells[:index])) / (nsites * len(taxa))
		utils.info("trimmed genes:", index)
		utils.info("trimmed sites:", nsites)
		utils.info("trimmed occupancy: %.2f" % cum_occupancy[index-1])
		utils.info("trimmed cell occupancy: %.2f" % cell_occupancy)

	diagnostics.log("nsites", nsites)
	diagnostics.log("cell_occupancy", cell_occupancy)

	# Write out supersequence for each taxon

	fasta = os.path.join(outdir, 'supermatrix.fa')
	with open(fasta, 'w') as f:
		for taxon in taxa:
			print >>f, ">%s" % taxon
			print >>f, ''.join(superseq[taxon][:index])

	partition = os.path.join(outdir, 'supermatrix.partition.txt')
	with open(partition, 'w') as f:
		cum_index = 0
		for i in xrange(0, index):
			print >>f, '%s, %s = %d-%d' % (
				partition_prefix, genes[i].name, cum_index+1, cum_index+sites[i])
			cum_index += sites[i]

	occupancy_txt = os.path.join(outdir, 'supermatrix.occupancy.txt')
	with open(occupancy_txt, 'w') as f:
		print >>f, '\t'.join(taxa)
		for row in occupancy[:index]:
			print >>f, '\t'.join(map(str, map(int, row)))

	cum_occupancy_txt = os.path.join(outdir, 'supermatrix.cum_occupancy.txt')
	with open(cum_occupancy_txt, 'w') as f:
		print >>f, '\n'.join(map(str, cum_occupancy[:index]))

	diagnostics.log('fasta', fasta)
	diagnostics.log('partition', partition)
	diagnostics.log('occupancy', occupancy_txt)
	diagnostics.log('cum_occupancy', cum_occupancy_txt)

# vim: noexpandtab ts=4 sw=4
