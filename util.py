import numpy as np
import networkx as nx
import re
import rsa
import scipy.stats

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
def vertex_disjoint_transform(G):
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
'''
def vertex_disjoint_paths(G,s,t,retrace=False):
	#TODO: figure out how many nodes are in the paths this algorithm would give and compare with current blacklisting etc algorithms
	#first modify G so that two-in two-out motifs evaluate correctly
	Gp = vertex_disjoint_transform(G)
	#then run max flow on that graph with the caveat that if we used the fork node transform on s we need to change the start to s
	if 'f{}'.format(s) in Gp.nodes:
		s = 'f{}'.format(s)

	R = nx.algorithms.flow.edmonds_karp(Gp,s,t)
	if retrace:
		paths = retrace_max_flow_paths(R,s,t)
		assert(len(paths) == R.graph['flow_value'])
		return paths
	else:
		return R.graph['flow_value']


'''
Given a residual network output of some nx flow algorithm, give me a list of <min vertex cut> VD paths 
'''
def retrace_max_flow_paths(R:nx.DiGraph,s,t):
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
			#keep moving

		#add on all neighbors *with flow == 1*
		for neigh in R.neighbors(current):
			if (neigh not in seen) and (R[current][neigh]['flow'] > 0):
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
def retrace_single_path(pred,current,path=None):
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
def str_compress(x:str):
	current = x[0]
	count = 1
	r = x[0]
	for i in range(1,len(x)):
		if (current != x[i]) or (i == len(x) - 1):
			if count > 3:
				r += '-' + (hex(count)[2:]) + '-' + x[i]
			else:
				r += r[-1] * (count - 1) + x[i]
			#the dashes and intermediate number add 3 chars total, so it only helps (when compared to just repeating the digit) if the count is more than 4 (it's the same if count = 4)
			current = x[i]
			count = 1
		else:
			count += 1

	return r