import numpy as np
import networkx as nx
import re
import rsa
import scipy.stats
from typing import List, Union

'''
'add one' (LE) to the set (thinking of it is a bitvector/binary set)

universe is a list to maintain ordering
'''
def set_increment(universe:list,current:set):
	for elt in universe:
		if elt not in current:
			current.add(elt)
			return current
		else:
			current.remove(elt)

	return current

def manh_norm(v):
	return sum(map(lambda x: abs(x),v))

def eucl_norm(v):
	return np.linalg.norm(v)


'''
transform G st vanilla max flow will count the paths correctly

this we do by adding, for each node with in-degree > 1 and out-degree >1, a dummy node (called a "fork" node)
that the original node points to and moving all of the original node's out-edges to be the out-edges of the fork node
'''
def vertex_disjoint_transform(G:nx.DiGraph) -> nx.DiGraph:
	Gp = nx.DiGraph()
	Gp.add_edges_from(G.edges)#copy G

	for v in list(Gp.nodes):#cast to list so that there are no concurrent modification issues
		#for every vertex
		if (Gp.in_degree(v) > 1) and (Gp.out_degree(v) > 1):#with in and out degree > 1 (forking vertex)
			#do the motif transform
			fork_node = 'f{}'.format(v)
			#assign all of v's out-edges to f* and remove them from v
			for u in list(Gp.neighbors(v)):#cast to list so that there are no concurrent modification issues
				Gp.add_edge(fork_node,u)
				Gp.remove_edge(v,u)
			#point v to the fork node
			Gp.add_edge(v,fork_node)

	#assign capacity 1 to all edges
	nx.set_edge_attributes(Gp,1,name='capacity')

	return Gp

'''
use max flow on a modified version of G to find the number of vertex disjoint paths between s and t

between the transform and the actual max flow algorithm, on big graphs this can take a looooooooooong time to run
'''
def vertex_disjoint_paths(G:nx.Graph,s,t,retrace=False) -> Union[List,int]:
	#first modify G so that two-in two-out motifs evaluate correctly
	Gp = vertex_disjoint_transform(nx.DiGraph(G))#shallow copy is fine
	#then run max flow on that graph with the caveat that if we used the fork node transform on s we need to change the start to s
	if 'f{}'.format(s) in Gp.nodes:
		s = 'f{}'.format(s)

	R = nx.algorithms.flow.edmonds_karp(Gp,s,t)
	if retrace:
		paths = retrace_max_flow_paths(R,s,t)
		assert(len(paths) == R.graph['flow_value'])
		return paths
	else:
		return int(R.graph['flow_value'])


'''
Given a residual network output of some nx flow algorithm, give me a list of <min vertex cut> VD paths 
'''
def retrace_max_flow_paths(R:nx.DiGraph,s,t) -> List:
	paths = []
	exp_q = [s]#exploration queue
	current = None
	seen = set()
	in_path = set()
	pred = {s:None}#map nodes to predecessors (this should be able to change)

	while len(exp_q) > 0:
		#run a DFS along some path to t, once we reach t, add that path to paths
		current = exp_q.pop(-1)
		if current == t:
			#add the path to paths
			paths.append(retrace_single_path(pred,t))
			#reset the things
			in_path = in_path.union(paths[-1]) - {t}
			seen = set()
			exp_q = [s]
			#keep moving

		#add on all neighbors *with flow == 1*
		for neigh in R.neighbors(current):
			if (neigh not in seen) and (neigh not in in_path) and (R[current][neigh]['flow'] > 0):
				if neigh in pred:
					pred[neigh] = current
				else:
					pred.update({neigh:current})
				exp_q.append(neigh)
				seen.add(neigh)

	return paths

'''
traverse the predecessor map to figure out a single path
'''
def retrace_single_path(pred,current,path=None) -> list:
	if path is None:
		path = []
	if pred[current] is None:
		return path

	fork_match = re.match(r'f([0-9]+)',str(current))
	if fork_match is not None:#don't add fork-nodes to the path
		assert(pred[current] == int(fork_match.group(1)))
		return retrace_single_path(pred,pred[current],path)
	else:
		return retrace_single_path(pred,pred[current],[current] + path)

'''
Compresses a hex string representation of some integer down by indicating repeated digits (mostly useful for human-readability)
'''
def str_compress(x:str) -> str:
	current = x[0]
	count = 1
	r = x[0]
	for i in range(1,len(x)):
		if (current != x[i]) or (i == len(x) - 1):
			if count > 3:
				r += '-' + (hex(count)[2:]) + '-' + x[i]
			else:
				r += r[-1] * (count - 1) + x[i]
			#the dashes, intermediate number, and preceding number add 4 chars total, so it only helps (when compared to just repeating the digit) if the count is more than 4 (it's the same if count = 4)
			current = x[i]
			count = 1
		else:
			count += 1

	return r

'''
undoes the compression from str_compress (note: outputs a string of hex)
'''
def str_decompress(x:str) -> str:
	r = x[0]
	repeat = 0#how many times should I repeat the next number?
	get_repetition_number = False
	rnum = ''
	for i in range(1,len(x)):
		if get_repetition_number:
			if x[i] == '-':
				#then we've gotten the whole thing
				repeat = int(rnum,16)-1
				rnum = ''
				get_repetition_number = False
			else:
				rnum += x[i]
		elif x[i] == '-':
			get_repetition_number = True
		else:
			r += (r[-1] * repeat) + x[i]
			repeat = 0

	if rnum != '':
		repeat = int(rnum,16)-1

	return r + (r[-1] * repeat)

'''
networkx-ized version of the older proprietary methods

uses the social network random generation algorithm (connected variant of https://www.pnas.org/content/99/suppl_1/2566)
'''
def generate_connected_rand_graph_from_deg_dist(num_nodes:int,approx_reciprocity=1.,distrib=lambda:scipy.stats.truncnorm.rvs(0,float('inf'),loc=3,scale=3)) -> nx.DiGraph:
	G = nx.DiGraph()
	G.add_node(0)

	#each node gets a degree from the distribution that is at least 1
	degrees = [max(1,int(distrib())) for _ in range(num_nodes)]
	connections_left_to_make = {i for i in range(num_nodes)}#these nodes all have connections left to make
	total_connections_left_to_make = sum(degrees)

	#randomly connect the nodes already in the network according to their degrees
	connections_made = 0
	while (connections_made < total_connections_left_to_make) and (len(connections_left_to_make) > 1):
		#pick a random starter node among the nodes already in the network
		i = np.random.choice(list(connections_left_to_make.intersection(G.nodes)))
		#pick a random other node to connect to i
		assign_to = np.random.choice(list(connections_left_to_make - {i} - set(G.neighbors(i))))
		if assign_to not in G.nodes:
			G.add_node(assign_to)

		G.add_edge(i,assign_to)
		connections_made += 1
		degrees[i] -= 1
		if degrees[i] == 0:
			connections_left_to_make.remove(i)

		#reciprocity "guarantees" from a given number
		if np.random.uniform(0,1) < approx_reciprocity:
			#then add the edge in the other direction
			G.add_edge(assign_to,i)
			connections_made += 1
			degrees[assign_to] -= 1
			if degrees[assign_to] == 0:
				connections_left_to_make.remove(assign_to)

	return G


'''
Generate a connected, undirected erdos-renyi random graph. Connect the graph by fudging vertices that aren't in the largest CC
'''
def generate_connected_ER_graph(num_nodes:int,avg_degree_approx:float,seed=None) -> nx.Graph:
	approx_p = avg_degree_approx / (2*(num_nodes - 1))#from that kh = 2*|E|/|V| and Ex[|E|] = n*(n-1)*p
	G = nx.generators.fast_gnp_random_graph(num_nodes,approx_p,seed=seed)
	np.random.seed(seed)
	#find the largest connected component (or rather, which nodes are not in it)
	ccs = [set(x) for x in sorted(nx.connected_components(G),key=len,reverse=False)]#biggest last

	#randomly connect using as few edges as possible all of these other ccs
	list_largest_cc = list(ccs[-1])
	for cc in ccs[:-1]:
		#pick a random node to connect to the main cc
		node_i = np.random.choice(list(cc))
		#pick a random node to connect this node to within the main cc
		node_j = np.random.choice(list_largest_cc)
		#connect them -- this cc is now a part of the main one
		G.add_edge(node_i,node_j)

	return G


'''
Calculate a (networkx) graph which is the union of all of these paths (which are just lists of things [anything hashable will do]). Assumes edges are undirected
'''
def path_union(paths:list) -> nx.Graph:
	G = nx.Graph()

	for path in paths:
		prev_node = path[0]
		for node in path[1:]:
			G.add_edge(prev_node,node)
			prev_node = node

	return G

'''
Given a list of lists of things (list of paths), find the maximum VD paths through the graph defined by those paths
'''
def vd_paths_from_candidates(paths:list,s,t) -> list:
	#first compute the path union
	U = path_union(paths)
	#then run the maxflow based (centralized) VD algorithm on that
	paths = vertex_disjoint_paths(U,s,t,retrace=True)
	return paths