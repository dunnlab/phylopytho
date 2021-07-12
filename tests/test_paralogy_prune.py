from phylopytho import treeprune
from dendropy import Tree
# import dendropy

def parse_tree( tree_string ):
	# Convert string to tree object, abstract this since it changes between dendropy versions
	return Tree.get_from_string( tree_string, schema='newick', rooting='force-rooted' )

def get_taxon_names( tree ):
	# Get frozenset of taxon names in a Tree
	return frozenset(treeprune.get_taxa(tree.leaf_nodes()))

def get_taxon_sets( trees ):
	# Take a list of trees
	taxa_sets = set()
	for tree in trees:
		taxa_sets.add( get_taxon_names( tree ) )
	return taxa_sets

def prune_to_sets( tree ):
	# paralogy_prune a tree, return a set of frozensets of taxon names in subtrees
	pruned_trees = []
	treeprune.paralogy_prune(tree, pruned_trees)
	return get_taxon_sets( pruned_trees )

def string_to_pruned_sets( tree_string ):
	# paralogy_prune a newick tree strong, return a set of frozensets of taxon names in subtrees
	return prune_to_sets( parse_tree( tree_string ) )

def test_tree1():
	# Two subtrees with 4 and 3 taxa

	# two reciprocally monophyletic orthologs
	tree1a = parse_tree ( "((A@1,(B@1,(C@1,D@1))), (A@2,(B@2,C@2)));" )
	tree1a_sets = prune_to_sets( tree1a )

	# same as above, different rooting
	tree1b = parse_tree ( "((B@1,(C@1,D@1)), (A@1,(A@2,(B@2,C@2))));" )
	tree1b_sets = prune_to_sets( tree1b )

	# same as above, different rooting
	tree1c = parse_tree ( "((C@1,D@1), (B@1,(A@1,(A@2,(B@2,C@2)))));" )
	tree1c_sets = prune_to_sets( tree1c )

	tree1_expected = set([frozenset(['A@2', 'C@2', 'B@2']), frozenset(['B@1', 'C@1', 'A@1', 'D@1'])])
	assert( tree1a_sets == tree1_expected )
	assert( tree1b_sets == tree1_expected )
	assert( tree1c_sets == tree1_expected )


def test_tree2():
	# Check that if the two maximally inclusive groups are the same size and overlap, they are pruned into three subtrees
	# ie, the intersection is its own tree

	tree2 = parse_tree ( "(A@1,(B@1,(C@1,(D@1,(A@2,B@2)))));" )
	tree2_sets = prune_to_sets( tree2 )
	tree2_expected = set([frozenset(['A@1', 'B@1']), frozenset(['C@1', 'D@1']), frozenset(['A@2', 'B@2'])])
	assert( tree2_sets == tree2_expected )

def test_tree3():
	# Build on test above by adding sister taxon that makes max groups different sizes, and make sure that the root is traversed when counting max groups
	tree3 = parse_tree ( "(E@1,(A@1,(B@1,(C@1,(D@1,(A@2,B@2))))));" )
	tree3_sets = prune_to_sets( tree3 )
	tree3_expected = set([frozenset(['B@1', 'C@1', 'A@1', 'D@1', 'E@1']), frozenset(['A@2', 'B@2'])])
	assert( tree3_sets == tree3_expected )

def test_tree4():
	# An example from Jesus Ballesteros Chavez
	# His original tree text: ((((((Sp1|A,Sp2|A),Sp3|A),Sp1|B),Sp2|B),Sp3|B),Sp1|C));
	tree4 = parse_tree ( "((((((A@1,B@1),C@1),A@2),B@2),C@2),A@3);" )
	tree4_sets = prune_to_sets( tree4 )

	tree4_expected = set([frozenset(['C@2', 'B@2', 'A@3']), frozenset(['B@1', 'C@1', 'A@1']), frozenset(['A@2'])])
	assert( tree4_sets == tree4_expected )

def test_tree5():
	# As above, but with an extra A so that there aren't two adjacent edges to prune
	tree5 = parse_tree ( "(((((((A@1,B@1),C@1),A@2),A@4),B@2),C@2),A@3);" )
	tree5_sets = prune_to_sets( tree5 )

	tree5_expected = set([frozenset(['C@2', 'B@2', 'A@3']), frozenset(['B@1', 'C@1', 'A@1']), frozenset(['A@2']), frozenset(['A@4'])])
	assert( tree5_sets == tree5_expected )

def test_tree6():
	# As above, but with another extra A
	tree6 = parse_tree("((((((((A@1,B@1),C@1),A@2),A@4),A@5),B@2),C@2),A@3);" )
	tree6_sets = prune_to_sets( tree6 )
	
	tree6_expected = set([frozenset(['A@1', 'B@1', 'C@1']), frozenset(['C@2', 'B@2', 'A@3']), frozenset(['A@2']), frozenset(['A@4']), frozenset(['A@5'])])
	assert( tree6_sets == tree6_expected )

def test_tree7():
	# As above, but with another extra A
	tree7 = parse_tree("(((((((((A@1,B@1),C@1),A@2),A@4),A@5),A@6),B@2),C@2),A@3);" )
	tree7_sets = prune_to_sets( tree7 )

	tree7_expected = set([frozenset(['A@1', 'B@1', 'C@1']), frozenset(['C@2', 'B@2', 'A@3']),frozenset(['A@2']), frozenset(['A@4']), frozenset(['A@5']), frozenset(['A@6'])])
	assert( tree7_sets == tree7_expected )

def test_tree8():
	# A tree with only sequences from the same species
	tree8 = parse_tree("(((((A@1,A@2),A@3),A@4),A@5),A@6);" )
	tree8_sets = prune_to_sets( tree8 )

	tree8_expected = set([frozenset(['A@1']), frozenset(['A@2']),frozenset(['A@3']), frozenset(['A@4']), frozenset(['A@5']), frozenset(['A@6'])])
	assert( tree8_sets == tree8_expected )

def test_tree9():
	# A tree with sequences all from different species
	tree9 = parse_tree("(((((A@1,B@1),C@1),D@1),E@1),F@1);" )
	tree9_sets = prune_to_sets( tree9 )

	tree9_expected = set([frozenset(['A@1', 'B@1', 'C@1', 'D@1', 'E@1', 'F@1'])])
	assert( tree9_sets == tree9_expected )


# Tests of slit_tree()
def test_split1():
	tree = parse_tree("(((((A@1,B@1),C@1),D@1),E@1),F@1);" )
	node1 = tree.mrca(taxon_labels=[ 'A@1','B@1' ])
	subtrees = treeprune.split_tree( tree, [ node1 ] )
	tree_sets = get_taxon_sets( subtrees )
	expected = set([frozenset([ 'A@1','B@1' ]), frozenset(['C@1', 'D@1', 'E@1', 'F@1'])])
	assert( tree_sets == expected )



def test_split2():
	tree = parse_tree("(((((A@1,B@1),C@1),D@1),E@1),F@1);" )
	node1 = tree.mrca(taxon_labels=[ 'A@1','B@1' ])
	node2 = tree.find_node_with_taxon_label('C@1')
	subtrees = treeprune.split_tree( tree, [ node1, node2 ] )
	tree_sets = get_taxon_sets( subtrees )
	expected = set([frozenset([ 'A@1','B@1' ]), frozenset(['C@1']), frozenset(['D@1', 'E@1', 'F@1'])])
	assert( tree_sets == expected )

def test_split3():
	tree = parse_tree("(((((A@1,B@1),C@1),D@1),E@1),F@1);" )
	node1 = tree.mrca(taxon_labels=[ 'A@1','B@1' ])
	node2 = tree.find_node_with_taxon_label('C@1')
	node3 = tree.mrca(taxon_labels=[ 'A@1','E@1' ])
	subtrees = treeprune.split_tree( tree, [ node1, node2, node3 ] )
	tree_sets = get_taxon_sets( subtrees )
	expected = set([frozenset([ 'A@1','B@1' ]), frozenset(['C@1']), frozenset(['D@1', 'E@1']), frozenset(['F@1'])])
	assert( tree_sets == expected )

def test_split4():
	tree = parse_tree("(((((A@1,B@1),C@1),D@1),E@1),F@1);" )
	node1 = tree.mrca(taxon_labels=[ 'A@1','B@1' ])
	node2 = tree.find_node_with_taxon_label('C@1')
	node3 = tree.mrca(taxon_labels=[ 'A@1','E@1' ])
	node4 = tree.find_node_with_taxon_label('F@1')
	subtrees = treeprune.split_tree( tree, [ node1, node2, node3 ] )
	tree_sets = get_taxon_sets( subtrees )
	expected = set([frozenset([ 'A@1','B@1' ]), frozenset(['C@1']), frozenset(['D@1', 'E@1']), frozenset(['F@1'])])
	assert( tree_sets == expected )
