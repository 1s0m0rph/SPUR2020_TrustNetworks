"""
Framework for testing VD algorithms

#TODO ADD TESTING (UNIT TESTS)!!! (where?)

#### Software engineering stuff#TODO

Function that takes in:
	- graph
	- algorithm to use
And outputs:
	- all these measures of how well it did

Script will: generate some graphs, run some algorithms on them, and output some numbers or plots
Example plot: GraphType of size 10, 100, 1000  (horizontal axis)
           vs Alg1 #paths_found, Alg2, ...
another for amount of compute (total? per node?)
(by the way: error bars)
another for number of nodes that need to participate?

In particular: ideally, a function that generates these graphs given number of desired nodes (and edges?)
	- "Erdos-Renyi"
	- Power-law or "social network" style
	- Grid sort of graph
	- Graphs we throw in for fun
"""
from typing import List, Union
from TN import *
from Hyperbolic import *
from stepper_sim import *

'''
graph generation wrapper

takes in an algorithm (function) as well as arguments for that function (which should generally at least contain num_nodes) and the type of node to produce (bust be a subclass of TN)
algorithm should produce an nx graph which we will then turn into a TN graph
'''
def generate_random_graph(algorithm:function,node_type,*args,**kwargs):
	nxG = algorithm(args,kwargs)
	tnG = convert_nx_graph_to_TN(nxG,node_type,args,kwargs)
	return tnG


'''
convert a given networkx graph or digraph to a list of TNNodes (or subclasses) that are the same graph
'''
def convert_nx_graph_to_TN(nxg:Union[nx.Graph,nx.DiGraph],node_type,*args,**kwargs):
	retG = []
	idx_map = {node:i for i,node in enumerate(nxg.nodes)}

	#nodes only require node_id
	if node_type in {TNNode_Stepper,TNNode}:
		retG = [node_type(i) for i in range(len(nxg.nodes))]

		for i,node in enumerate(nxg.nodes):
			for neighbor in nxg.neighbors(node):
				retG[i].add_public_key_in_person(retG[idx_map[neighbor]])

	#nodes require additional parameters to construct
	elif node_type == HyperNode:
		retG = [node_type(i,args,kwargs) for i in range(len(nxg.nodes))]

		#hyperbolic graphs must be built in a specific order (treelike) -- use a BFS tree to construct
		root = nxg.nodes[0]
		bfst = nx.algorithms.traversal.bfs_tree(nxg,root)#bfs tree from networkx
		#use the tree to construct an ordering (pre-order tree traversal [parent then children])
		node_order = nx_bfs_tree_preordering(bfst,root)

		for i,node in enumerate(node_order):
			for neighbor in nxg.neighbors(node):
				retG[i].add_public_key_in_person(retG[idx_map[neighbor]])


	return retG

'''
construct a (pre) ordering of the nodes in bfst
'''
def nx_bfs_tree_preordering(bfst,current,order=None):
	if order is None:
		order = []

	#pre-order
	order = [current] + order
	#call on daughters
	for d in bfst.neighbors(current):
		order = nx_bfs_tree_preordering(bfst,d,order)

	return order


'''
Run vd_path_alg on G from s to t

path alg possibilities:
	'TN-v2':	 		TNNode.count_vd_paths_to_v2(destination id, return paths = True)
	'TN-v3':			TNNode.count_vd_paths_to_v3(destination id, destination coordinates, heuristic='dot')
	'synch-v2':			TNNode_Stepper.count_vd_paths_to_v2(destination id, return paths = True, TTL = infinity)
	'hyper-addr':		HyperNode.count_vd_paths_to_hyper_from_addr(dest address, npaths=inf, max distance scale = inf, stop on first fail = false)
	'hyper':			HyperNode.count_vd_paths_to_hyper(dest coordinates, npaths=inf, max distance scale = inf, stop on first fail = false)
	'hyper-multi-addr':	HyperNode.count_vd_paths_to_hyper_multibl_from_addr(dest address, max distance scale = inf, stop on first fail = false)
	'hyper-multi':		HyperNode.count_vd_paths_to_hyper_multibl(dest coordinates, max distance scale = inf, stop on first fail = false)
	'local-mf':			hyper_VD_paths_local(hyperbolic graph, start, target, max distance scale = inf, distance measure = 'path', autoscale increment = none)


returns:
	success: bool: true iff all paths found were valid and all paths found were VD. if validate=False, this is always true. if no paths are found, this is always true
	num_paths_found: how many paths did we find?
	num_nodes_used: how many nodes had to *send a message* -- nothing else counts towards this
	total_messages_sent: how many messages were sent?
	
does NOT return (i.e. these need to be calculated later):
	number of possible paths
	lower-bound network usage
	

this function will reset the graph after calculating the needed quantities
'''
def run_single_pair(G:List[TNNode],vd_path_alg:str,s:int,t:int,validate=True):#TODO
	pass


'''
As with run_single_pair, but this runs many and calculates summary stats for the metrics

if random_sample is set, we will randomly sample among the possible pairs of paths without replacement min(max_npairs,floor(num_possible_pairs*sample_proportion)) pairs to calculate
'''
def run_many_pairs(G:List[TNNode],vd_path_alg:str,random_sample=False,sample_proporiton=1.,max_npairs=float('inf'),validate=True):#TODO
	pass

