from unittest import TestCase
from TNUtil import *

def graphs_are_same(nxg:nx.Graph,tng:List[TNNode]):
	for node in nxg.nodes:
		for node_j in nxg.nodes:
			if nxg.has_edge(node,node_j):
				if tng[node_j] not in tng[node].neighbors:
					return False
			else:
				if tng[node_j] in tng[node].neighbors:
					return False

	return True

class TestExperimentRunner(TestCase):
	def test_convert_nx_graph_to_TN(self):
		nxg = nx.Graph()
		nxg.add_edges_from([[0,1],[0,2],[1,2],[2,3],[2,4],[1,4]])
		tng = convert_nx_graph_to_TN(nxg,TNNode)
		tng_nx = convert_to_nx_graph(tng,nx.Graph)
		assert(nx.is_isomorphic(nxg,tng_nx))

		#directed graph
		nxg = nx.DiGraph()
		nxg.add_edges_from([[0,1],[0,2],[1,2],[2,3],[2,4],[1,4]])
		tng = convert_nx_graph_to_TN(nxg,TNNode)
		tng_nx = convert_to_nx_graph(tng)
		assert(nx.is_isomorphic(nxg,tng_nx))

		#hypernodes
		nxg = nx.Graph()
		nxg.add_edges_from([[0,2],[0,3],[0,4],[1,3],[1,4],[1,5],[2,3],[3,4],[3,5],[3,6],[5,6]])
		tng = convert_nx_graph_to_TN(nxg,HyperNode,3)
		#make sure all the coords are set up
		for node in tng:
			assert((type(node.coords) == np.complex128) or (type(node.coords) == complex))
		tng_nx = convert_to_nx_graph(tng,nx.Graph)
		assert(nx.is_isomorphic(nxg,tng_nx))


	def test_nx_bfs_tree_preordering(self):
		#small graph for debugging
		# nxg = nx.Graph()
		# nxg.add_edges_from([[0,2],[0,3],[0,4],[1,3],[1,4],[1,5],[2,3],[3,4],[3,5],[3,6],[5,6]])
		#massive random graph for testing
		num_nodes = 100
		nxg = generate_connected_ER_graph(num_nodes,0.1)
		bfst = nx.algorithms.traversal.bfs_tree(nxg,0)
		node_order = nx_bfs_tree_preordering(bfst,0)

		#first node in the order should be 0
		assert(node_order[0] == 0)

		#single contiguous connected component test
		nodes_in_graph = {node_order.pop(0)}
		for node in node_order:
			edge_exists = False
			for ng in nodes_in_graph:
				if nxg.has_edge(node,ng):
					edge_exists = True
					break
			assert(edge_exists)
			nodes_in_graph.add(node)

		#all nodes should be present
		assert(nodes_in_graph == set(nxg.nodes))


	def test_generate_random_graph(self):
		#mostly want to make sure args work how we expect them to
		generate_random_graph(generate_connected_ER_graph,TNNode,[100,0.5],{'seed':0})
		generate_random_graph(generate_connected_ER_graph,HyperNode,[100,0.5],{'seed':0},[3])

	def test_pick_n_random_pairs(self):
		#small graph for debugging
		nxg = nx.Graph()
		nxg.add_edges_from([[0,2],[0,3],[0,4],[1,3],[1,4],[1,5],[2,3],[3,4],[3,5],[3,6],[5,6]])
		G =  convert_nx_graph_to_TN(nxg,HyperNode,3)#q=3
		#this graph has 10 possible pairs
		npairs = 10
		pairs = pick_n_random_pairs(G,npairs)
		assert(len(pairs) == npairs)

		for s,t in pairs:
			assert(G[t] not in G[s].neighbors)

		npairs_fail = 11
		self.assertRaises(Exception,pick_n_random_pairs,G,npairs_fail)