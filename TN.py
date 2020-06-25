"""
General utils etc for TN/Trust Networks
"""

from util import *

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

	def reset_flags(self):
		self.resetted_flag = True
		self.pulse_pred = {}
		self.search_blacklist_flag = False

	'''
	do another graph search to reset pulse numbers
	'''
	def reset_search(self):
		if self.resetted_flag:
			#we're done here
			return

		#otherwise set our predecessor, reset our number, and tell all our neighbors to reset
		self.reset_flags()

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

def generate_rand_graph_from_deg_dist(num_nodes,dist=lambda:scipy.stats.truncnorm.rvs(0,float('inf'),loc=3,scale=3),approx_reciprocity=1.,node_type=TNNode):#TODO remove
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


def convert_to_nx_graph(TNG,graph_type=nx.DiGraph):
	nxG = graph_type()
	nxG.add_nodes_from([node.id for node in TNG])
	for node in TNG:
		for connection in node.neighbors:
			nxG.add_edge(node.id,connection.id,capacity=1)

	return nxG