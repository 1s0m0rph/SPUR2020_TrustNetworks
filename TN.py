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
		self.operations_done = 0  #how many messages have I had to *send* during the current run?#TODO rename this to messages_sent or something

		self.address = None#and will stay that way unless annotated by an Embedder
		self.adist = None#this is just the distance function provided by an Embedder

	def __repr__(self):
		if self.address is not None:
			return 'Trust Network Node {} with address {}'.format(self.id,self.address)
		else:
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
		
		
	def reset_all_and_return_ops(self):
		self.reset_flags()
		ops = self.operations_done
		self.operations_done = 0
		return ops

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

		self.operations_done += 1#in order to retrace these paths, we need to send a message
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


	def vertex_blacklist_search_greedy(self,target_addr,max_paths=float('inf'),max_dist_scale=float('inf'),stop_on_first_failure=False):
		if self.address is None:
			raise TypeError("Trying to run a greedy search on a graph which has not been completely addressed")

		self.search_blacklist_flag = True
		st_dist = self.adist(self.address,target_addr)
		max_dist = st_dist*max_dist_scale

		#begin the search
		neighbors_to_call = sorted([(self.adist(n.address,target_addr),n) for n in self.neighbors],key=lambda x: x[0])
		paths = []
		pnum = 0
		for dist,neighbor in neighbors_to_call:
			if dist <= max_dist:
				self.operations_done += 1
				path_ret = neighbor.vertex_bl_greedy_interm(target_addr,self,max_dist,stop_on_first_failure,pnum)
				if path_ret is not None:
					paths.append(path_ret)
					if len(paths) > max_paths:
						break
				elif stop_on_first_failure:
					break
			pnum += 1

		return paths

	def vertex_bl_greedy_interm(self,target_addr,pred,max_dist,stop_on_first_failure,pulse_num):
		if self.address == target_addr:
			self.pulse_pred.update({pulse_num:pred})
			self.resetted_flag = False
			self.search_blacklist_flag = False
			path = self.v2_vd_blacklist_zip(pulse_num,[])
			return path

		if pulse_num in self.pulse_pred:
			return None

		if self.search_blacklist_flag:
			return None

		self.pulse_pred.update({pulse_num:pred})
		self.resetted_flag = False

		neighbors_to_call = sorted([(self.adist(n.address,target_addr),n) for n in self.neighbors],key=lambda x: x[0])
		for dist,neighbor in neighbors_to_call:
			if dist <= max_dist:
				self.operations_done += 1
				path_ret = neighbor.vertex_bl_greedy_interm(target_addr,self,max_dist,stop_on_first_failure,pulse_num)
				if path_ret is not None:
					return path_ret
				#TODO does it make sense to propagate this?
				elif stop_on_first_failure:#small algorithm change: propagate this
					return None

		return None

def convert_to_nx_graph(TNG,graph_type=nx.DiGraph):
	nxG = graph_type()
	nxG.add_nodes_from([node.id for node in TNG])
	for node in TNG:
		for connection in node.neighbors:
			nxG.add_edge(node.id,connection.id,capacity=1)

	return nxG