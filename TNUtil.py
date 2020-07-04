"""
Slightly more specialized than util.py (requires TNNode and subclasses), but not so specialized as to really belong in run_experiment.py
"""

from util import *
from Hyperbolic import *
from stepper_sim import *

'''
graph generation wrapper

takes in an algorithm (function) as well as arguments for that function (which should generally at least contain num_nodes) and the type of node to produce (bust be a subclass of TN)
algorithm should produce an nx graph which we will then turn into a TN graph

we have multiple args/kwargs pairs so unfortunately they need to be done all fudge-like like this
'''
def generate_random_graph(algorithm,node_type,algorithm_args=(),algorithm_kwargs=None,node_args=(),node_kwargs=None) -> List[TNNode]:
	if node_kwargs is None:
		node_kwargs = {}
	if algorithm_kwargs is None:
		algorithm_kwargs = {}
	nxG = algorithm(*algorithm_args,**algorithm_kwargs)
	tnG = convert_nx_graph_to_TN(nxG,node_type,*node_args,**node_kwargs)
	return tnG


'''
convert a given networkx graph or digraph to a list of TNNodes (or subclasses) that are the same graph
'''
def convert_nx_graph_to_TN(nxg:Union[nx.Graph,nx.DiGraph],node_type,*args,**kwargs) -> List[TNNode]:
	retG = []
	nodes = list(nxg.nodes)
	idx_map = {node:i for i,node in enumerate(nodes)}#maps node names from the original graph onto indexes in retG

	#nodes only require node_id
	if node_type in [TNNode_Stepper,TNNode]:
		retG = [node_type(i) for i in range(len(nodes))]

		for i,node in enumerate(nodes):
			for neighbor in nxg.neighbors(node):
				retG[i].add_public_key_in_person(retG[idx_map[neighbor]])

	#nodes require additional parameters to construct
	elif node_type == HyperNode:
		retG = [node_type(i,*args,**kwargs) for i in range(len(nodes))]

		#hyperbolic graphs must be built in a specific order (treelike) -- use a BFS tree to construct
		root = nodes[0]
		bfst = nx.algorithms.traversal.bfs_tree(nxg,root)#bfs tree from networkx
		#use the tree to construct an ordering (pre-order tree traversal [parent then children])
		node_order = nx_bfs_tree_preordering(bfst,root)
		retG[idx_map[node_order[0]]].init_as_root()

		for node in node_order:
			for neighbor in nxg.neighbors(node):
				retG[idx_map[node]].add_public_key_in_person(retG[idx_map[neighbor]])


	return retG

'''
construct a (pre) ordering of the nodes in bfst
'''
def nx_bfs_tree_preordering(bfst,current,order=None):
	if order is None:
		order = []

	#pre-order
	order = order + [current]
	#call on daughters
	for d in bfst.neighbors(current):
		order = nx_bfs_tree_preordering(bfst,d,order)

	return order

def pick_n_random_pairs(G:List[Union[TNNode,HyperNode,TNNode_Stepper]],npairs:int):
	#use numpy instead of trying to be too fancy
	pairs = [(s,t) for s in range(len(G)) for t in range(s,len(G)) if (G[t] not in G[s].neighbors) and (s != t)]
	return [pairs[i] for i in np.random.choice(list(range(len(pairs))),size=npairs,replace=False)]