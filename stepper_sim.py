"""
General idea: we want to figure out how well our algorithms work on an actual network with distributed parallel computation. We could literally build such a network, but that would be difficult and likely time-consuming. Thus, instead, we'll build a time-step based simulation to run computations concurrently in a sort of simulation
"""

from TN import *

#it's rewind time

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

class TNNode_Stepper(TNNode):

	"""
	A version of TNNode that is compatible with stepper-simulations

	general function format:
		needs to return a list of all the function calls generated by the function (fn,args,kwargs) as well as a callback for processing once all subcalls have returned
	"""

	def __init__(self,node_id):
		super().__init__(node_id)

		self.time = 0#this will be used to make sure no function calls happen prematurely
		self.operations = []#queue of operations we need to perform (fn,args,kwargs)
		self.paths = []#maintained at the object level to make callbacks easier
		self.neighbors_to_call = []
		self.pulse_num = 0

	def __repr__(self):
		return 'SSTNN {}'.format(self.id)

	'''
	move our self-stored time variable up by 1 and do all of the operations in our queue
	'''
	def increment_time(self):
		# self.time += 1
		this_step_operations = self.operations.copy()
		self.operations.clear()
		for fn, args, kwargs, callback in this_step_operations:
			ret = fn(*args, **kwargs)
			callback(ret)


	def count_vd_paths_to_v2_callback(self,path):
		if path is not None:
			self.paths.append(path)
			if len(self.neighbors_to_call) > 0:
				self.pulse_num += 1
				self.neighbors_to_call.pop(0).v2_vd_paths_interm_synchronized(path[-1].id,self,self.pulse_num,self.time + 1,self.count_vd_paths_to_v2_callback)

	def cleanup(self):
		# now reset pulse numbers and blacklist flags
		for n in self.neighbors:
			n.reset_search()

		self.search_blacklist_flag = False

	'''
	2nd naive version of vertex-disjoint path count algorithm between s (this/self) and t
	this version runs until no paths are found and can return the found paths for verification of VD property
	'''
	def count_vd_paths_to_v2(self,dest_id,return_paths=True):
		self.search_blacklist_flag = True
		#start a search to dest from each neighbor
		pulse_num = 0
		self.neighbors_to_call = list(self.neighbors)
		self.neighbors_to_call.pop(0).v2_vd_paths_interm_synchronized(dest_id,self,pulse_num,self.time+1,self.count_vd_paths_to_v2_callback)


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
	def v2_vd_paths_interm_synchronized(self,dest_id,pred,pulse_num,time,callback):
		if time > self.time:
			#then we need to postpone processing this part until our time is updated
			self.operations.append((self.v2_vd_paths_interm_synchronized,[dest_id,pred,pulse_num,time,callback],{},callback))
			return None

		if pulse_num in self.pulse_pred:
			return None#we've already been visited in this search

		if self.id == dest_id:
			#blacklist nodes on this path
			self.pulse_pred.update({pulse_num:pred})
			path = self.v2_vd_blacklist_zip(pulse_num,[])
			return path  #for now we're just counting paths, but we could also reconstruct what they are from this point

		if self.search_blacklist_flag:
			return None  #we are already blacklisted (or have been visited in this search), so don't go this way

		#we've relayed this pulse now
		self.pulse_pred.update({pulse_num:pred})

		#otherwise ask all our neighbors if they know the muffin man
		for n in self.neighbors:
			path = n.v2_vd_paths_interm_synchronized(dest_id,self,pulse_num,time+1,callback)
			if path is not None:
				return path

		return None

	'''
	I am on the path to dest from s, so I need to be blacklisted

	this method also reconstructs the path, mostly for verification
	'''
	def v2_vd_blacklist_zip(self,pulse_num,path: list):
		if pulse_num not in self.pulse_pred:
			return path  #this could also return the reconstructed path, but would require more data to be passed between nondes

		self.search_blacklist_flag = True
		return self.pulse_pred[pulse_num].v2_vd_blacklist_zip(pulse_num,[self] + path)

#how many paths does this method find on average, as a percentage of the maximum (empirical)?

N = 100

np.random.seed(0)
tng = generate_rand_graph_from_deg_dist(N,approx_reciprocity=1,node_type=TNNode_Stepper)

prop_sum = 0
npairs = 0
pair_count = 0

for s in range(N):
	for t in range(s+1,N):
		if t not in tng[s].neighbors:
			npairs += 1#doing this first so we can get progress reports

# s = 0
# t = 8

for s in range(N):
	for t in range(s+1,N):
		if t not in tng[s].neighbors:
			if pair_count % 100 == 0:
				print('{} of {} ({:.3f}%)'.format(pair_count,npairs,100.*float(pair_count)/float(npairs)))
			exact_total_paths = vertex_disjoint_paths(convert_to_nx_graph(tng),s,t)

			tng[s].count_vd_paths_to_v2(t)
			all_done_flag = False
			while not all_done_flag:
				all_done_flag = True
				for tn in tng:
					tn.time += 1  # increment first so all of the nodes are incremented before we do this time step
				for tn in tng:
					tn.increment_time()
				for tn in tng:
					if len(tn.operations) > 0:  # has to be separate because high node ids can give operations to low ones
						all_done_flag = False

			tng[s].cleanup()  # this would probably just be done either by a final pulse from s or by timeout

			paths = tng[s].paths

			# print('{} of {} total paths found: '.format(len(paths),exact_total_paths))

			nodes_seen = set()

			for path in paths:
				prev_node = tng[s]
				for node in path:
					# verify that paths exist
					if node not in prev_node.neighbors:
						print('EDGE DOES NOT EXIST: ({},{})'.format(prev_node.id,node.id))
					prev_node = node
					# verify that paths are disjoint
					if (node.id != s) and (node.id != t) and (node.id in nodes_seen):
						print('REPEATED NODE ID: {}'.format(node.id))
					else:
						nodes_seen.add(node.id)
				# print(path)

			prop_sum += float(len(paths)) / float(exact_total_paths)
			pair_count += 1

			#reset all the graph things (this could be done in reality with either a pulse or a timeout)
			for i in range(len(tng)):
				tng[i].pulse_pred = {}
				tng[i].resetted_flag = False
				tng[i].search_blacklist_flag = False
				tng[i].time = 0
				tng[i].operations = []
				tng[i].paths = []
				tng[i].neighbors_to_call = []
				tng[i].pulse_num = 0

print('AVERAGE PROPORTION OF PATHS FOUND: {:.3f}%'.format(100*(prop_sum/npairs)))