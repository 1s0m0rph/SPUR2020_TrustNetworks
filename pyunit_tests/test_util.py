from unittest import TestCase
from util import *


class TestGenerate_connected_ER_graph(TestCase):
	def test_generate_connected_ER_graph(self):
		#first check that a high-p graph is just connected end of story
		a = generate_connected_ER_graph(100,0.7,seed=0)
		assert(nx.is_connected(a))#regression test

		#now check that it successfully connects a graph that is almost surely not connected
		a = generate_connected_ER_graph(100,0.01,seed=0)
		assert(nx.is_connected(a))

	def test_vertex_disjoint_transform(self):
		G = nx.DiGraph()
		G.add_edges_from([[0,1],[0,2],[1,3],[2,3],[3,4],[3,5],[4,6],[5,6]])#double-diamond graph (2-in, 2-out at node 3)
		Gp = vertex_disjoint_transform(G)
		Gp_exp = nx.DiGraph()
		Gp_exp.add_edges_from([[0,1],[0,2],[1,3],[2,3],[3,'f3'],['f3',4],['f3',5],[4,6],[5,6]])

		assert(Gp.edges == Gp_exp.edges)

	def test_vertex_disjoint_paths(self):
		G = nx.DiGraph()
		G.add_edges_from([[0,1],[0,2],[1,3],[2,3],[3,4],[3,5],[4,6],[5,6]])  # double-diamond graph (2-in, 2-out at node 3)
		#number of paths only
		npaths = vertex_disjoint_paths(G,0,6,retrace=False)
		assert(npaths == 1)

		#large random graph for retracing (tests subfunctions as well as this one)
		G = generate_connected_ER_graph(100,0.3,0)
		#random start, random (non-neighbor) end
		s = np.random.choice(list(G.nodes))
		t = np.random.choice(list(set(G.nodes) - set(G.neighbors(s))))
		#get paths
		paths = vertex_disjoint_paths(G,s,t,retrace=True)
		#check that paths are valid (i.e. no repeated nodes, all edges exist, and goes from s to t) and VD
		nodes_seen = set()
		for path in paths:
			prev_node = s#starts at s
			for node in path[:-1]:
				assert(node not in nodes_seen)#paths are VD (+ no repeated nodes/edges)
				nodes_seen.add(node)
				assert(G.has_edge(prev_node,node))#all edges exist
				prev_node = node
			assert(G.has_edge(path[-2],path[-1]))#final edge exists
			assert(path[-1] == t)#ends at t

	def test_str_compress(self):
		x = 'EB000000000041D8'
		comp = str_compress(x)
		assert(comp == 'EB0-a-41D8')

		x = 'DD00044EEEEE0012'
		comp = str_compress(x)
		assert(comp == 'DD00044E-5-0012')

		x = '011111000111111A'
		comp = str_compress(x)
		assert(comp == '01-5-0001-6-A')

	def test_str_decompress(self):
		x = 'EB0-a-41D8'
		comp = str_decompress(x)
		assert (comp == 'EB000000000041D8')

		x = 'DD00044E-5-0012'
		comp = str_decompress(x)
		assert (comp == 'DD00044EEEEE0012')

		x = '01-5-0001-6-A'
		comp = str_decompress(x)
		assert (comp == '011111000111111A')