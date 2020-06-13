"""
General utils etc for TN/Trust Networks
"""

import numpy as np
import networkx as nx
import rsa
import scipy.stats

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
def vertex_disjoint_paths(G,s,t):
	#first modify G so that two-in two-out motifs evaluate correctly
	Gp = vertex_disjoint_transform(G)
	#then run max flow on that graph with the caveat that if we used the fork node transform on s we need to change the start to s
	if 'f{}'.format(s) in Gp.nodes:
		s = 'f{}'.format(s)

	return nx.algorithms.flow.edmonds_karp(Gp,s,t).graph['flow_value']

class TNNode:
	INITIAL_MAX_NEIGHBOR_SEPARATION = 1#this is k, in the math
	MAX_NEIGHBOR_SEPARATION_INCREMENT = 1#for dynamic separation increase, what should the increment size be?
	NUM_COORD_DIMENSIONS = 2
	ADDRESS_BIT_SIZE = 256

	def __init__(self,node_id):
		self.id = node_id
		self.private_key = None
		self.public_key = None
		self.neighbors = set()
		self.neighbors_ids = set()

		self.is_eve = False

		self.coords = np.zeros((self.NUM_COORD_DIMENSIONS,),dtype=np.uint16)#this will be set when we first add a neighbor
		#FIXME coordinate vector dtypes need to be dynamic with number of coordinate dimensions

		#self.generate_keys()#this would be done in reality, but we won't do it now for the sake of determinism

		#algorithmic data
		#are we blacklisted for the current search (in a real implementation this would need to be a list)
		self.search_blacklist_flag = False
		#for the current pulse (in real system: in a given search), who was our predecessor for each given pulse?
		self.pulse_pred = {}#maps pulse nums onto pred objects
		#visited flag used for resetting pulse numbers
		self.resetted_flag = False#in a real system this would be a search id that would automatically reset after some time for each node

		# metrics
		self.operations_done = 0  # have I had to do any work in the current pass through?

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

		neighbor_coord_centroid = np.mean(np.array([n.coords for n in self.neighbors],dtype=np.uint16),axis=0,dtype=np.uint16)
		max_separation = self.INITIAL_MAX_NEIGHBOR_SEPARATION
		disp,choices = self.generate_coords(max_separation,self.NUM_COORD_DIMENSIONS,-np.array(neighbor_coord_centroid,dtype=np.int32))
		#check if these coords are okay by pinging the network
		while self.check_is_coord_in_use_init(neighbor_coord_centroid+disp):
			choices_left = True
			for coord_choices in choices:
				if len(coord_choices) == 0:#all have to have valid choices, not just some
					choices_left = False

			if not choices_left:
				#increase k
				max_separation += self.MAX_NEIGHBOR_SEPARATION_INCREMENT
				choices = None#for the generation method

			disp,choices = self.generate_coords(max_separation,self.NUM_COORD_DIMENSIONS,-np.array(neighbor_coord_centroid,dtype=np.int32),choices)


		self.coords = neighbor_coord_centroid + disp

	'''
	generate n_dim coords such that the norm (manhattan so it's guaranteed to be integers) of the vector formed by those coords is leq max_norm

	coord choices is used to ensure determinism
	max negative is used to ensure the coordinates are always non-negative
		this is a list of length n_dim with integer value <= 0 for each element
	'''
	def generate_coords(self,max_norm,n_dim,max_negative:list,coord_choices=None):
		single_bound = int(max_norm / n_dim)#this is manhattan-norm specific; norm of biggest vector will be this times n_dim which will be no larger than max_norm
		if coord_choices is None:
			coord_choices = [{i for i in range(max(-single_bound,max_negative[j]),single_bound+1)} for j in range(n_dim)]#want to be able to go negative
		coords = []
		for i in range(n_dim):
			coord = np.random.choice(list(coord_choices[i]))
			#don't make this choice again
			coord_choices[i].remove(coord)
			coords.append(coord)

		return np.array(coords,dtype=np.int32),coord_choices

	def check_is_coord_in_use_init(self,coords):
		self.pulse_pred.update({-1:None})
		self.resetted_flag = False
		is_in_use = False
		for neighbor in self.neighbors:
			if neighbor.check_is_coord_in_use(coords,self,-1):
				is_in_use = True
				break

		self.reset_search()
		return is_in_use

	'''
	semicentralized (fully decentralized is a trivial step from here) check for in-use coordinates (so coords are unique)
	'''
	def check_is_coord_in_use(self,coords,pred,pulse_num):
		if pulse_num in self.pulse_pred:
			return False

		if np.allclose(self.coords,coords):
			return True

		self.pulse_pred.update({pulse_num:pred})
		self.resetted_flag = False#this is for later

		#otherwise ask our neighbors if they're using that number
		for neighbor in self.neighbors:
			in_use = neighbor.check_is_coord_in_use(coords,self,pulse_num)
			if in_use:
				return True

		return False



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
			return 1#I have 1 path to myself (just a base case so everything works out fine)

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
			return path#for now we're just counting paths, but we could also reconstruct what they are from this point

		if self.search_blacklist_flag or (pulse_num in self.pulse_pred):
			return None#we are already blacklisted (or have been visited in this search), so don't go this way

		#we've relayed this pulse now
		self.pulse_pred.update({pulse_num:pred})
		self.resetted_flag = False

		#otherwise ask all our neighbors if they know the muffin man
		for n in self.neighbors:
			path = n.v2_vd_paths_interm(dest_id,self,pulse_num,)
			if path is not None:
				return path

		return None

	'''
	I am on the path to dest from s, so I need to be blacklisted

	this method also reconstructs the path, mostly for verification
	'''
	def v2_vd_blacklist_zip(self,pulse_num,path: list):
		if pulse_num not in self.pulse_pred:
			return path#this could also return the reconstructed path, but would require more data to be passed between nondes

		self.search_blacklist_flag = True
		return self.pulse_pred[pulse_num].v2_vd_blacklist_zip(pulse_num,[self] + path)


	"""
	V3
	
	(doesn't need to be synchronized since it's depth first and therefore won't be done in parallel)
	"""

	'''
	clean up the predecessor tree
	'''
	def pred_tree_clean(self):
		#unset predecessor (preprocess)
		self.pulse_pred = {}

		successors = []
		for neighbor in self.neighbors:
			if (-1 in neighbor.pulse_pred) and (neighbor.pulse_pred[-1] == self):
				successors.append(neighbor)

		#call pred tree clean on all successors
		for succ in successors:
			succ.pred_tree_clean()

	'''
	3rd version uses heuristics -- coordinates
	
	if you want the heuristic function to be sorted backwards (that is, in descending order), negate the output of the heuristic function
	
	built-in heuristics:
		'dot' (default) : largest among neighbors n dot product of v-n vector and v-t vector ('most correct direction')
		'dist' : shortest distance to t among all neighbors n ('closest neighbor') -- min over n of norm(coord(n) - coord(t)), using manhattan norm
		'dist-manh': synonym for dist
		'dist-eucl': as with dist-manh, but uses euclidean norm
		
		so far all the heuristics seem to perform basically the same (except that manhattan distance is faster)
	'''
	def count_vd_paths_to_v3(self,dest_id,dest_coords,heuristic=None):
		if (heuristic is None) or (heuristic == 'dot'):
			heuristic = lambda x: - np.dot(self.coords - x.coords, self.coords - dest_coords)#heuristic A: largest dot-product (most in the right direction)
		elif (heuristic == 'dist') or (heuristic == 'dist-manh'):
			heuristic = lambda x: manh_norm(x.coords - dest_coords)
		elif heuristic == 'dist-eucl':
			heuristic = lambda x:eucl_norm(x.coords - dest_coords)

		self.search_blacklist_flag = True
		#start a search to dest from each neighbor
		neighbors_to_call = list(sorted(list(self.neighbors),key=heuristic))
		paths = []
		for neighbor in neighbors_to_call:
			self.operations_done += 1
			path = neighbor.v3_vd_paths_interm(dest_id,dest_coords,self,heuristic)
			if path is not None:
				paths.append(path)
				self.pred_tree_clean()#TODO does it make sense to pred clean here/is there somewhere else it should be done?

		for neighbor in self.neighbors:
			neighbor.reset_search()

		self.search_blacklist_flag = False

		return paths

	'''
	someone has asked me to find dest
	'''
	def v3_vd_paths_interm(self,dest_id,dest_coords,pred,heuristic):
		self.operations_done += 1
		if self.id == dest_id:#this has to happen before the visited check to implement path shadowing
			#blacklist nodes on this path
			self.pulse_pred.update({-1:pred})
			self.resetted_flag = False
			self.search_blacklist_flag = False#t is never blacklisted, but the function after sets it to be so
			path = self.v2_vd_blacklist_zip(-1,[])
			return path

		if -1 in self.pulse_pred:
			return None#we've already been visited in this search

		if self.search_blacklist_flag:
			return None#we are already blacklisted (or have been visited in this search), so don't go this way

		#we've relayed this pulse now
		self.pulse_pred.update({-1:pred})
		self.resetted_flag = False

		#otherwise ask the *right* neighbor(s) if they know the muffin man
		neighbors_to_call = list(sorted(list(self.neighbors),key=heuristic))
		for neighbor in neighbors_to_call:
			path = neighbor.v3_vd_paths_interm(dest_id,dest_coords,self,heuristic)
			if path is not None:
				return path

		return None


def generate_rand_graph_from_deg_dist(num_nodes,dist=lambda:scipy.stats.truncnorm.rvs(0,float('inf'),loc=3,scale=3),approx_reciprocity=1.,node_type=TNNode):
	#use the social network rand gen algorithm to make a new network
	G = [node_type(i) for i in range(num_nodes)]

	#assign a degree to each node
	degrees = [max(1,int(dist())) for i in range(num_nodes)]#make sure everything is connected to at least one
	connections_left_to_make = {i for i in range(num_nodes)}#these nodes all have connections left to make
	total_connections_to_make = sum(degrees)

	#randomly connect the nodes according to their degrees
	connections_made = 0
	while (connections_made < total_connections_to_make) and (len(connections_left_to_make) > 1):
		#pick a random starter node
		i = np.random.choice(list(connections_left_to_make))
		#pick a random other node to assign to
		assign_to = G[np.random.choice(list(connections_left_to_make - {i} - G[i].neighbors))]#no self-edges, no multiedges
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
