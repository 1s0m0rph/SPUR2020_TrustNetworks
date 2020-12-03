import numpy as np
import random
import networkx as nx
import re
import scipy.stats
from typing import List, Union, Tuple
# import sys

# sys.setrecursionlimit(10000)#big graphs require this

PROGRESS_INTERVAL = 0

def is_iterable(t):
	try:
		x = iter(t)
	except TypeError:
		return False
	return True

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

cache the zero-flow residual network so nx doesn't have to calculate it every time (this results in a 2.1x speedup!)
'''
def vertex_disjoint_paths(G:nx.Graph,s,t,retrace=False,Gp=None) -> Union[List,int]:
	#first modify G so that two-in two-out motifs evaluate correctly
	if Gp is None:
		Gp = nx.algorithms.flow.build_residual_network(vertex_disjoint_transform(nx.DiGraph(G)),'capacity')#shallow copy is fine
	#else reset (doesn't need to be done since the networkx functions all do that for us)
	#then run max flow on that graph with the caveat that if we used the fork node transform on s we need to change the start to s
	if 'f{}'.format(s) in Gp.nodes:
		s = 'f{}'.format(s)

	R = nx.algorithms.flow.preflow_push(Gp,s,t,residual=Gp)#TODO maybe more optimizations here -- still looking at 67% of the time being spent in this fn
	if retrace:
		paths = retrace_max_flow_paths(R,s,t)
		return paths
	else:
		return int(R.graph['flow_value'])


'''
Given a residual network output of some nx flow algorithm, give me a list of <min vertex cut> VD paths 
'''
def retrace_max_flow_paths(R:nx.DiGraph,s,t) -> List:
	paths = [[]]
	exp_q = [s]#exploration queue
	current = None
	seen = set()
	in_path = set()

	while len(paths) < R.graph['flow_value']+1:#this is the correct number
		#run a DFS along some path to t, once we reach t, add that path to paths
		current = exp_q.pop(-1)
		if current == t:
			#reset the things
			in_path = in_path.union(paths[-1]) - {t}
			seen = set()
			exp_q = [s]
			paths.append([])
			#keep moving

		#add on all neighbors *with flow == 1*
		for neigh in R.neighbors(current):
			if (neigh not in seen) and (neigh not in in_path) and (R[current][neigh]['flow'] > 0):
				if (type(neigh) != str) or (neigh[0] != 'f'):#fork nodes should not be in paths since they only exist in the augmented network
					paths[-1].append(neigh)#no need for predecessors since we're being greedy
				exp_q.append(neigh)
				seen.add(neigh)
				break#small optimization

	return paths[:-1]#we always add a dummy at the end so drop it

'''
traverse the predecessor map to figure out a single path
'''
def retrace_single_path(pred,current,path=None) -> list:#TODO remove deprecated
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


"""
GRAPH GENERATORS
"""


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
	approx_p = avg_degree_approx / (num_nodes - 1)#from that kh = (2*|E|)/|V| and Ex[|E|] = (n*(n-1)/2)*p
	G = nx.generators.fast_gnp_random_graph(num_nodes,approx_p,seed=seed)
	np.random.seed(seed)
	#find the largest connected component (or rather, which nodes are not in it)
	ccs = [list(x) for x in sorted(nx.connected_components(G),key=len,reverse=False)]#biggest last

	#randomly connect using as few edges as possible all of these other ccs
	list_largest_cc = ccs[-1]
	for cc in ccs[:-1]:
		#pick a random node to connect to the main cc
		node_i = np.random.choice(cc)
		#pick a random node to connect this node to within the main cc
		node_j = np.random.choice(list_largest_cc)
		#connect them -- this cc is now a part of the main one
		G.add_edge(node_i,node_j)

	return G


'''
variable-dense graph: 
	kh = s*f(|V|) for constant scale s, kh = Theta(f(|V|), |E| = Theta(|V| f(|V|))

e.g. if f(x) = log x, then |E| = Theta(|V| log |V|)

function types:
	<any function from float to float>
	'log': ln x
	'sqrt': sqrt x
	'dense': x (standard dense graph)
	'sparse': 1 (standard sparse graph)
'''
def generate_connected_variable_dense_ER_graph(num_nodes:int,scale:float,ftype:str,seed=None) -> nx.Graph:
	f = None
	if ftype == 'log':
		f = np.log
	elif ftype == 'sqrt':
		f = np.sqrt
	elif ftype == 'dense':
		f = lambda x: x
	elif ftype == 'sparse':
		f = lambda x: 1.
	else:
		raise AttributeError("Unknown density function type: {}".format(ftype))

	avg_deg = scale * f(num_nodes)
	return generate_connected_ER_graph(num_nodes,avg_deg,seed=seed)

def generate_nbad_unioning_fast(num_nodes:int) -> nx.Graph:
	# first find the nearest n
	n = int(np.round((np.sqrt(4 * num_nodes - 11) + 1) / 2))
	# initialize the graph
	G = nx.Graph()
	G.add_nodes_from(['s','t'])
	#start by making the zeroth (s-t) path
	prev_node = 's'
	for node_idx in range(n-1):
		G.add_edge(prev_node,(0,node_idx))
		prev_node = (0,node_idx)
	G.add_edge(prev_node,'t')

	path_max_intersected = [(2*i)-1 for i in range(n-1)]#what is the maximum index this path has had to intersect with another at?

	for path_num in range(1,n):
		#make the ith path
		num_nodes_in_path = n + path_num - (1 if path_num != (n-1) else 0) #not necessarily the number of nodes that will be added
		prev_node = 's'
		#make the dummy nodes
		num_dummies = path_num
		for i in range(num_dummies):
			this_node = (path_num,i)
			G.add_edge(prev_node,this_node)
			prev_node = this_node

		#connect up all of our intersections starting from path 0, node i-1 (total intersections = i)
		intersections_to_make = path_num
		for with_path in range(intersections_to_make):
			#intersect this path at node given by the max_intersected array
			path_max_intersected[with_path] += 1
			intersect_with_at = (with_path,path_max_intersected[with_path])
			G.add_edge(prev_node,intersect_with_at)
			prev_node = intersect_with_at

		#straight on til morning
		next_node = (path_num,num_dummies+intersections_to_make)
		for _ in range(num_nodes_in_path - (num_dummies + intersections_to_make)):
			G.add_edge(prev_node,next_node)
			prev_node = next_node
			next_node = (path_num,next_node[1]+1)

		#final edge goes to t
		G.add_edge(prev_node,'t')

	return G

'''
Number of nodes is only approximate -- the real number is the nearest num_nodes st num_nodes = n^2 - n + 3 for some integer n

guaranteed to have ceil(n/2) vd paths between s and t

interesting note: these graphs appear to always be planar
'''
def generate_nbad_unioning(num_nodes:int) -> nx.Graph:
	#first find the nearest n
	n = int(np.round((np.sqrt(4*num_nodes - 11) + 1)/2))
	#initialize the graph
	G = nx.Graph()
	G.add_nodes_from(['s','t'])
	#now construct all of the paths in a rudimentary way between s and t
	for i in range(n):
		prev_node = 's'
		for j in range(n + i - 1):
			G.add_edge(prev_node,(i,j))
			prev_node = (i,j)
		if i != n-1:
			G.add_edge(prev_node,'t')
		else:
			G.add_edges_from([(prev_node,(i,n+i-2)),((i,n+i-2),(i,n+i-1)),((i,n+i-1),'t')])#dummy node at the end

	#now construct the merge table
	M = [[() for _ in range(n + i - 1)] for i in range(n)]
	#now fill in the merge table
	#start with base cases
	for i in range(n-1):
		M[0][i] = (i+1,i+1)
		M[i+1][i+1] = (0,i)
		#merge these
		G = nx.contracted_nodes(G,(0,i),(i+1,i+1),self_loops=False)
	ip,jp = (2,3)
	for i in range(1,n-1):
		for j in range(i+1,n+i-1):
			if len(M[i][j]) == 0:
				M[i][j] = (ip,jp)
				M[ip][jp] = (i,j)
				#do the merge now
				G = nx.contracted_nodes(G,(i,j),(ip,jp),self_loops=False)
				ip += 1
				jp += 1
		ip = i + 2
		jp = 2*(i+2) - 1

	#remove self-loops since for some reason contracted nodes doesn't always fail to add them
	for node in G:
		if G.has_edge(node,node):
			G.remove_edge(node,node)

	return G

'''
n-bad example for neighbor-blacklist, assuming the pathfinding algorithm finds shortest paths

number of nodes as a function of n: 5 + 4(n-1) [(s,t,s1,m,t1,t) + (si,pi,l1i,l2i,l3i,ti for i in 1,n)]
so n as function of num_nodes = ((num_nodes - 5)/6) + 1
'''
def generate_nbad_neighborbl(num_nodes:int) -> nx.Graph:
	n = int(np.round(((num_nodes-5)/4)+1))
	G = nx.Graph()

	#add s and t and all s_i,t_i; connect s to all s_i and t to all t_i
	G.add_nodes_from(['s','t'])
	for i in range(n):
		G.add_edges_from([('s','s_{}'.format(i)),('t','t_{}'.format(i))])

	#finish path 1 (shortest)
	G.add_edges_from([('s_0','m'),('m','t_0')])

	#add all the other paths
	for i in range(1,n):
		si = 's_{}'.format(i)
		ti = 't_{}'.format(i)
		l1i = 'l1_{}'.format(i)
		l2i = 'l2_{}'.format(i)
		G.add_edges_from([(si,l1i),(l1i,l2i),(l2i,ti),#path i [length 5]
						  (si,'m'),#path st_i [length 4]
						  ('m',ti)#the real killer
						  ])

	return G

'''
n-bad example for both neighbor-blacklist AND naive blacklisting/best-of-n-paths

number of nodes as a function of n: n^2 + 3n + 2
'''
def generate_nbad_either(num_nodes:int):
	n = int(np.round((np.sqrt(4*num_nodes + 1) - 3)/2))

	nnodes_U = n**2 - n + 3
	U = generate_nbad_unioning_fast(nnodes_U)

	nnodes_M = 4*(n-1) + 5
	M = generate_nbad_neighborbl(nnodes_M)

	M.remove_node('t')
	U.remove_node('s')

	G = nx.union(M,U,rename=('M-','U-'))

	#merge the appropriate nodes
	for i in range(n):
		#connect 'M-t_{i}' with 'U-({i},0)'
		G.add_edge('M-t_{}'.format(i),'U-({}, 0)'.format(i))

	G = nx.relabel_nodes(G,{'M-s':'s','U-t':'t'})
	return G

'''
this function assumes it's always an edgelist (maybe add others if needed)
first arg just here for the sake of regularity
'''
def read_graph_from_data(_,filename) -> nx.Graph:
	G = nx.read_edgelist(filename,create_using=nx.Graph,delimiter=',',nodetype=int)
	assert(nx.is_connected(G))#if not this won't work
	return G

def type_check(obj,t):
	if is_iterable(obj) and is_iterable(t):
		for et,at in zip(t,obj):
			if is_iterable(et):
				if type(at) not in et:
					raise TypeError("Expected address type {}, got {}".format(et,at))
			else:
				if et != type(at):
					raise TypeError("Expected address type {}, got {}".format(et,at))
	elif type(obj) != t:
		if is_iterable(t):
			if type(obj) not in t:
				raise TypeError("Expected address type {}, got {}".format(t,obj))
		else:
			raise TypeError("Expected address type {}, got {}".format(t, type(obj)))