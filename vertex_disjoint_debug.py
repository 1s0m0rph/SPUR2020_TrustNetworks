from TN import *

'''
transform G st vanilla max flow will count the paths correctly

this we do by adding, for each node with in-degree > 1 and out-degree >1, a dummy node (called a "fork" node)
that the original node points to and moving all of the original node's out-edges to be the out-edges of the fork node
'''
def vertex_disjoint_transform(G):
	Gp = nx.DiGraph()
	Gp.add_edges_from(G.edges)  #copy G

	for v in list(Gp.nodes):  #cast to list so that there are no concurrent modification issues
		#for every vertex
		if (Gp.in_degree(v) > 1) and (Gp.out_degree(v) > 1):  #with in and out degree > 1 (forking vertex)
			#do the motif transform
			fork_node = 'f{}'.format(v)
			#assign all of v's out-edges to f* and remove them from v
			for u in list(Gp.neighbors(v)):  #cast to list so that there are no concurrent modification issues
				Gp.add_edge(fork_node,u)
				Gp.remove_edge(v,u)
			#point v to the fork node
			Gp.add_edge(v,fork_node)

	#assign capacity 1 to all edges
	nx.set_edge_attributes(Gp,1,name='capacity')

	return Gp


'''
use max flow on a modified version of G to find the number of vertex disjoint paths between s and t
'''
def vertex_disjoint_paths(G,s,t):
	#first modify G so that two-in two-out motifs evaluate correctly
	Gp = vertex_disjoint_transform(G)
	#then run max flow on that graph with the caveat that if we used the fork node transform on s we need to change the start to s
	if 'f{}'.format(s) in Gp.nodes:
		s = 'f{}'.format(s)

	return nx.algorithms.flow.edmonds_karp(Gp,s,t).graph['flow_value']

s = 0
t = 5

#apparently networkx has an edmonds-karp implementation, let's try that
np.random.seed(0)
tng = generate_rand_graph_from_deg_dist(10,approx_reciprocity=0.5)

exact_total_paths = vertex_disjoint_paths(convert_to_nx_graph(tng),s,t)

paths = tng[s].count_vd_paths_to_v2(t)

print('{} of {} total paths found: '.format(len(paths),exact_total_paths))

for path in paths:
	print(path)