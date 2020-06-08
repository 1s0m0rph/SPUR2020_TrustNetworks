import numpy as np
import networkx as nx
import rsa
import scipy.stats

class TNNode:

	def __init__(self,node_id):
		self.id = node_id
		self.private_key = None
		self.public_key = None
		self.neighbors = set()
		self.neighbors_ids = set()

		self.is_eve = False

		#self.generate_keys()#this would be done in reality, but we won't do it now for the sake of determinism

		#algorithmic data
		#are we blacklisted for the current search (in a real implementation this would need to be a list)
		self.search_blacklist_flag = False
		#for the current pulse (in real system: in a given search), who was our predecessor for each given pulse?
		self.pulse_pred = {}#maps pulse nums onto pred objects
		#visited flag used for resetting pulse numbers
		self.resetted_flag = False#in a real system this would be a search id that would automatically reset after some time for each node

	def generate_keys(self):
		self.public_key,self.private_key = rsa.newkeys(16)

	def __repr__(self):
		return 'Trust Network Node with ID {}'.format(self.id)

	'''
	simulates high-fidelity in-person key transfer
	no trust assumptions yet
	'''
	def add_public_key_in_person(self,t_node):
		self.neighbors.add(t_node)
		self.neighbors_ids.add(t_node.id)

	'''
	tell me how many *paths* (not necessarily vertex-disjoint!) there are from me to some destination node
	naive implementation (assumes no eves)
	'''
	def count_paths_to(self,dest_id):
		#since this is a surrogate for trust, and we're just trying to find keys, infinite trust for already-known keys
		if self.id == dest_id:
			return float('inf')
		if dest_id in self.neighbors_ids:
			return float('inf')

		#now we just have to ask our neighbors how many paths there are
		paths = 0
		for neighbor in self.neighbors:
			paths += neighbor.path_count_interm(set(),dest_id)

		return paths

	'''
	I have been asked by a neighbor to count the number of paths to this person
	'''
	def path_count_interm(self,nodes_visited: set,dest_id):
		if dest_id == self.id:
			return 1  #I have 1 path to myself (just a base case so everything works out fine)

		#otherwise I'll have to ask my neighbors
		nodes_visited.add(self)
		paths = 0
		for neighbor in (self.neighbors - nodes_visited):
			paths += neighbor.path_count_interm(nodes_visited,dest_id)

		return paths

	'''
	2nd naive version of vertex-disjoint path count algorithm between s (this/self) and t
	this version runs until no paths are found and can return the found paths for verification of VD property
	'''
	def count_vd_paths_to_v2(self,dest_id,return_paths=True):
		self.search_blacklist_flag = True
		#start a search to dest from each neighbor
		path_count = 0
		pulse_num = 0
		if return_paths:
			paths = []
		for n in self.neighbors:
			path = n.v2_vd_paths_interm(dest_id,self,pulse_num)
			if path is not None:
				path_count += 1
				if return_paths:
					paths.append(path)
			pulse_num += 1

		#now reset pulse numbers and blacklist flags
		for n in self.neighbors:
			n.reset_search()

		self.search_blacklist_flag = False

		if return_paths:
			return [list(reversed(path)) for path in paths]#because of how they're reconstructed, the paths must be reversed first
		else:
			return path_count

	'''
	do another graph search to reset pulse numbers
	'''
	def reset_search(self):
		if self.resetted_flag:
			#we're done here
			return

		#otherwise set our predecessor, reset our number, and tell all our neighbors to reset
		self.resetted_flag = True
		self.pulse_pred = {}
		self.search_blacklist_flag = False

		for n in self.neighbors:
			n.reset_search()

	'''
	someone has asked me to find dest
	'''
	def v2_vd_paths_interm(self,dest_id,pred,pulse_num):
		if self.id == dest_id:
			#blacklist nodes on this path
			self.pulse_pred.update({pulse_num:pred})
			path = self.v2_vd_blacklist_zip(pulse_num,[])
			return path  #for now we're just counting paths, but we could also reconstruct what they are from this point

		if self.search_blacklist_flag or (pulse_num in self.pulse_pred):
			return None  #we are already blacklisted (or have been visited in this search), so don't go this way

		#we've relayed this pulse now
		self.pulse_pred.update({pulse_num:pred})

		#otherwise ask all our neighbors if they know the muffin man
		for n in self.neighbors:
			path = n.v2_vd_paths_interm(dest_id,self,pulse_num)
			if path is not None:
				return path

		return None

	'''
	I am on the path to dest from s, so I need to be blacklisted
	
	this method also reconstructs the path, mostly for verification
	'''
	def v2_vd_blacklist_zip(self,pulse_num,path:list):
		if pulse_num not in self.pulse_pred:
			return path#this could also return the reconstructed path, but would require more data to be passed between nondes

		self.search_blacklist_flag = True
		return self.pulse_pred[pulse_num].v2_vd_blacklist_zip(pulse_num,path + [self])


def generate_rand_graph_from_deg_dist(num_nodes,dist=lambda:scipy.stats.truncnorm.rvs(0,float('inf'),loc=3,scale=3),approx_reciprocity=1.):
	#use the social network rand gen algorithm to make a new network
	G = [TNNode(i) for i in range(num_nodes)]

	#assign a degree to each node
	degrees = [max(1,int(dist())) for i in range(num_nodes)]  #make sure everything is connected to at least one
	connections_left_to_make = {i for i in range(num_nodes)}  #these nodes all have connections left to make
	total_connections_to_make = sum(degrees)

	#randomly connect the nodes according to their degrees
	connections_made = 0
	while (connections_made < total_connections_to_make) and (len(connections_left_to_make) > 1):
		#pick a random starter node
		i = np.random.choice(list(connections_left_to_make))
		#pick a random other node to assign to
		assign_to = G[np.random.choice(list(connections_left_to_make - {i} - G[i].neighbors))]  #no self-edges, no multiedges
		G[i].add_public_key_in_person(assign_to)
		connections_made += 1
		degrees[i] -= 1
		if degrees[i] == 0:
			connections_left_to_make.remove(i)

		#reciprocity guarantees based on the given number
		if np.random.uniform(0,1) < approx_reciprocity:
			#then add the edge going back to i
			assign_to.add_public_key_in_person(G[i])
			connections_made += 1
			degrees[assign_to.id] -= 1
			if degrees[assign_to.id] == 0:
				connections_left_to_make.remove(assign_to.id)

	return G


def convert_to_nx_graph(TNG):
	nxG = nx.DiGraph()
	nxG.add_nodes_from([node.id for node in TNG])
	for node in TNG:
		for connection in node.neighbors:
			nxG.add_edge(node.id,connection.id,capacity=1)

	return nxG


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