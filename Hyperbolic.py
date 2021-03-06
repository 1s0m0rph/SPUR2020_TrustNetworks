from TN import *

"""
####Algorithm ideas

Orthogonal question: when to stop trying paths
- when you can't find any more? Seems to use too many nodes
- when there is a single failure? Seems to work, but still uses too many nodes
- maybe s tries each neighbor once, then stops.



Optimal performance looks not-so-good at this point. We may be able to fix it a bit using a (semifudgy) solution:
	Each node maintains a local view of the network out to k hops
	When a node receives a request to route a message to t, it runs some kind of more well-informed algorithm (how do local max-flows translate to s-t paths?), then sends the packet (along with its recommendations?) to some neighbor
	

dense networks seem to work better (as long as we're not looking for all the edges) what about semidense networks (kh = log n or sqrt(n) or something)?

Packets sent between s and t represent keyshares. Should we add an edge between s and t if they successfully communicate? Under what extra conditions (e.g. number of agreeing paths)?
"""

'''
arguments are complex numbers
'''
def hyper_dist(a:complex,b:complex):
	return np.arccosh(1 + (2*(abs(a-b)**2))/((1 - abs(a)**2)*(1 - abs(b)**2)))#FIXME why is this taking so long? Maybe we should precalculate/cache all neighbor distances

class Isometry:

	def __init__(self,rotation,translation):
		self.rot = rotation
		self.trans = translation

	def __repr__(self):
		return 'ISO r = {}; t = {}'.format(self.rot,self.trans)

	'''
	arg 0 is an isometry (pair of complex numbers), arg 1 is a single complex number
	'''

	def evaluate(self,arg):
		return (self.rot * arg + self.trans) / (1 + np.conj(self.trans) * self.rot * arg)

	def cross(self,l):
		return Isometry((self.rot * l.rot + l.rot * self.trans * np.conj(l.trans)) / (self.rot * l.trans * np.conj(self.trans) + 1),
						(self.rot * l.trans + self.trans) / (self.rot * l.trans * np.conj(self.trans) + 1))

	def inv(self):
		a = Isometry(np.conj(self.rot),0)
		b = Isometry(1,-self.trans)
		return a.cross(b)


'''
algorithm 1
'''
def define_generators(q):
	generators = []
	rot_isom = Isometry(np.e ** (1j * (2 * np.pi / q)),0)
	trans_isom = Isometry(1,np.tanh(np.arccosh(1 / (np.sin(np.pi / q))))).cross(Isometry(-1,0))

	for i in range(q):
		#for some reason doing it the way their pseudocode says to doesn't work because of this R^i thing
		#it only affects the zeroth generator (and therefore only the root)
		rot_isom_i = Isometry((np.array([rot_isom.rot,rot_isom.trans]) ** complex(i))[0],0j)
		generators.append(rot_isom_i.cross(trans_isom).cross(rot_isom_i.inv()))

	return generators

class HyperNode(TNNode):

	def __init__(self,node_id,q,coordinates=None,index=None,isometry=None):
		super().__init__(node_id)
		self.coords = coordinates
		self.idx = index
		self.isom = isometry
		self.q = q

		self.d_coords = []#list of available daughter coords (list of [(coords,index,isometry)] which can instantiate daughters)

		self.neighbors = set()
		self.daughters = set()

		self.distances = {}#map (cache) certain nodes to distance to that node
		#this causes hyper_dist to be called 40% as often -- ~1.1x speedup overall (impact of this function on the entire run is now halved)

		self.d_add_search_flag = False

		self.max_neighbor_called = -1#in a real system, one per search

		self.times_visited = 0

	def __repr__(self):
		return "HN, {} (@{})".format(self.id,self.coords)

	def reset_all_and_return_ops(self):
		ops = super().reset_all_and_return_ops()
		self.d_add_search_flag = False
		self.max_neighbor_called = -1
		return ops

	'''
	algorithm 2
	'''
	def init_as_root(self):
		self.coords = 0 + 0j
		self.idx = 0
		self.isom = Isometry(1 + 0j,0 + 0j)
		self.calculate_daughter_coords()

	'''
	algorithm 3
	'''
	def calculate_daughter_coords(self):
		generators = define_generators(self.q)

		for i in range(self.q):
			didx = (self.idx + i) % self.q
			disom = self.isom.cross(generators[didx])
			dcoord = disom.evaluate(0)

			self.d_coords.append((dcoord,didx,disom))

	def reset_flags(self):
		super().reset_flags()
		self.d_add_search_flag = False
		self.max_neighbor_called = -1
		self.distances.clear()#just so the memory doesn't fill up. this is probably no bueno but since the target node changes all the time there's no use in remembering previous targets

	'''
	add this node as a daughter and give it the info it needs (coords, index, isometry)
	'''
	def add_daughter(self,d,visited=None):
		if len(self.d_coords) == 0:
			self.calculate_daughter_coords()  #only do this if/when we add a daughter
		if visited is None:
			visited = set()
		if self in visited:
			return None
		visited.add(self)
		if len(self.d_coords) > 0:
			self.daughters.add(d)
			return self.d_coords.pop(0)
		else:
			#ask our neighbors if they can add d
			for neighbor in self.neighbors:
				info = neighbor.add_daughter(d,visited)
				if info is not None:
					return info
			return None

	'''
	add a link from me to n
	'''
	def add_neighbor(self,n):
		if n in self.neighbors:
			return#already exists
		if self.coords is None:
			if n.coords is None:
				raise AttributeError("Tried to connect two nodes that weren't already in the network")
			#make n our parent
			self.coords, self.idx, self.isom = n.add_daughter(self)

		self.neighbors.add(n)
		n.add_neighbor(self)#enforce reciprocity

	def add_public_key_in_person(self,t_node):
		self.add_neighbor(t_node)


	'''
	calculate the distance to this other node in the embedding space
	'''
	def dist_to(self,target_coords:complex):
		if target_coords in self.distances:
			return self.distances[target_coords]
		else:
			dist = hyper_dist(self.coords,target_coords)
			self.distances.update({target_coords:dist})#this shouldn't get too big since we're only calculating distance to a small set of target nodes
			return dist


	"""
	algorithm 4
	
	this is NOT as stated in the paper, since we're using blacklisting (is it guaranteed to still work? it certainly seems to)
	"""

	'''
	initialize the search
	
	npaths variable tells us how many paths to find (we'll stop when we find this many or when we have found the max).
	max distance scale tells us how far nodes are allowed to be from t, as a linear function of the distance between s and t
		specifically, nodes that are further than max_dist_scale * (dist from s to t) are excluded
	'''
	def count_vd_paths_to_hyper(self,dest_coords,max_paths=float('inf'),max_dist_scale=float('inf'),stop_on_first_failure=False):
		self.search_blacklist_flag = True
		st_dist = hyper_dist(self.coords,dest_coords)
		#start a search to dest from each neighbor
		neighbors_to_call = list(sorted(list(self.neighbors),key=lambda x: x.dist_to(dest_coords)))
		paths = []
		for neighbor in neighbors_to_call:
			if hyper_dist(neighbor.coords,dest_coords) <= max_dist_scale * st_dist:
				self.operations_done += 1
				path_ret = neighbor.gnh_interm(dest_coords,self,st_dist,max_dist_scale)
				if path_ret is not None:
					paths.append(path_ret)
					if len(paths) >= max_paths:
						break
				elif stop_on_first_failure:
					break

		return paths

	def gnh_interm(self,dest_coords,pred,st_dist,max_dist_scale):
		if self.coords == dest_coords:#this has to happen before the visited check to implement path shadowing
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
		neighbors_to_call = list(sorted(list(self.neighbors),key=lambda x: x.dist_to(dest_coords)))
		for neighbor in neighbors_to_call:
			if hyper_dist(neighbor.coords,dest_coords) <= max_dist_scale * st_dist:
				self.operations_done += 1
				path_ret = neighbor.gnh_interm(dest_coords,self,st_dist,max_dist_scale)
				if path_ret is not None:
					return path_ret

		return None
	
	
	"""
	neighbor-blacklisting
	"""

	'''
	initialize the search

	npaths variable tells us how many paths to find (we'll stop when we find this many or when we have found the max).
	max distance scale tells us how far nodes are allowed to be from t, as a linear function of the distance between s and t
		specifically, nodes that are further than max_dist_scale * (dist from s to t) are excluded
	'''
	def count_vd_paths_to_hyper_neighborbl(self,dest_coords,max_dist_scale=float('inf'),stop_on_first_failure=False):
		st_dist = hyper_dist(self.coords,dest_coords)
		#start a search to dest from ourselves
		candidate_paths = []
		pulse_num = 0
		while self.max_neighbor_called < (len(self.neighbors)-1):#semantically identical to just doing the neighbor loop
			path_ret = self.gnh_neighborbl_interm(dest_coords,None,pulse_num,st_dist,max_dist_scale)
			pulse_num += 1
			if path_ret is not None:
				candidate_paths.append([self] + path_ret)
			elif stop_on_first_failure:
				break

		#now calculate a graph union among all the candidate paths and use that to do a strict max flow
		if len(candidate_paths) > 0:
			paths = vd_paths_from_candidates(candidate_paths,self,candidate_paths[0][-1])
		else:
			paths = []

		return paths


	def gnh_neighborbl_interm(self,dest_coords,pred,pulse_num,st_dist,max_dist_scale):
		if self.coords == dest_coords:#this has to happen before the visited check to implement path shadowing
			#blacklist nodes on this path
			self.pulse_pred.update({pulse_num:pred})
			self.resetted_flag = False
			path = self.neighbor_blacklist_zip(pulse_num)
			return path

		if pulse_num in self.pulse_pred:
			return None#we've already been visited in this search

		#we've relayed this pulse now
		self.pulse_pred.update({pulse_num:pred})
		self.resetted_flag = False

		#otherwise ask the *right* neighbor(s) if they know the muffin man
		neighbors_to_call = list(sorted(list(self.neighbors),key=lambda x:x.dist_to(dest_coords)))
		start_idx = self.max_neighbor_called + 1
		for nidx in range(start_idx,len(neighbors_to_call)):
			neighbor = neighbors_to_call[nidx]
			if hyper_dist(neighbor.coords,dest_coords) <= max_dist_scale * st_dist:
				self.operations_done += 1
				path_ret = neighbor.gnh_neighborbl_interm(dest_coords,self,pulse_num,st_dist,max_dist_scale)
				self.max_neighbor_called = nidx
				if path_ret is not None:
					return path_ret

		self.max_neighbor_called = len(self.neighbors) - 1

		return None

	'''
	reconstruct the path from s to t backwards and blacklist the nodes on it
	
	iterative now to avoid recursion depth limits
	'''
	def neighbor_blacklist_zip(self,pulse_num):
		path = []
		current = self
		while (pulse_num not in current.pulse_pred) or (current.pulse_pred[pulse_num] is None):
			path = [current] + path
			current = current.pulse_pred[pulse_num]

		return path


	"""
	multi-visit algorithm
	
	generally, each node asks its neighbors how many times they've been visited and prioritizes them based on the ones which have been visited the fewest times
	"""

	def get_num_times_visited(self):
		self.operations_done += 1#this is why this method is important
		return self.times_visited

	'''
	initialize the search

	npaths variable tells us how many paths to find (we'll stop when we find this many or when we have found the max).
	max distance scale tells us how far nodes are allowed to be from t, as a linear function of the distance between s and t
		specifically, nodes that are further than max_dist_scale * (dist from s to t) are excluded
	'''
	def count_vd_paths_to_hyper_multivisit(self,dest_coords,max_dist_scale=float('inf'),stop_on_first_failure=False):
		st_dist = hyper_dist(self.coords,dest_coords)
		#start a search to dest from ourselves
		candidate_paths = []
		pulse_num = 0
		while self.max_neighbor_called < (len(self.neighbors) - 1):#semantically identical to just doing the neighbor loop
			path_ret = self.gnh_multivisit_interm(dest_coords,None,pulse_num,st_dist,max_dist_scale)
			pulse_num += 1
			if path_ret is not None:
				candidate_paths.append([self] + path_ret)
			elif stop_on_first_failure:
				break

		#now calculate a graph union among all the candidate paths and use that to do a strict max flow
		if len(candidate_paths) > 0:
			paths = vd_paths_from_candidates(candidate_paths,self,candidate_paths[0][-1])
		else:
			paths = []

		return paths

	def gnh_multivisit_interm(self,dest_coords,pred,pulse_num,st_dist,max_dist_scale):
		if self.coords == dest_coords:#this has to happen before the visited check to implement path shadowing
			#blacklist nodes on this path
			self.pulse_pred.update({pulse_num:pred})
			self.resetted_flag = False
			path = self.multi_visit_zip(pulse_num)
			return path

		if pulse_num in self.pulse_pred:
			return None#we've already been visited in this search

		#we've relayed this pulse now
		self.pulse_pred.update({pulse_num:pred})
		self.resetted_flag = False

		#otherwise ask the *right* neighbor(s) if they know the muffin man
		neighbors_blacklist_counts = {n:n.get_num_times_visited() for n in self.neighbors}
		self.operations_done += len(self.neighbors)#send all of the blacklist-count requests
		#primary sort is the blacklist count, secondary is distance to t
		neighbors_to_call = list(sorted(list(self.neighbors),key=lambda x:(neighbors_blacklist_counts[x],x.dist_to(dest_coords))))
		start_idx = self.max_neighbor_called + 1
		for nidx in range(start_idx,len(neighbors_to_call)):
			neighbor = neighbors_to_call[nidx]
			if hyper_dist(neighbor.coords,dest_coords) <= max_dist_scale * st_dist:
				self.operations_done += 1#send the pathfind request
				path_ret = neighbor.gnh_multivisit_interm(dest_coords,self,pulse_num,st_dist,max_dist_scale)
				self.max_neighbor_called = nidx
				if path_ret is not None:
					return path_ret

		self.max_neighbor_called = len(self.neighbors) - 1

		return None

	'''
	I am on the path to dest from s, so I need to be blacklisted

	this method also reconstructs the path
	
	changed to loop to avoid recursion depth problems
	'''
	def multi_visit_zip(self,pulse_num):
		#TODO implement tree-return for multi-blacklist zip?
		path = []
		current = self
		while (pulse_num in current.pulse_pred) and (current.pulse_pred[pulse_num] is not None):
			current.operations_done += 1#in order to retrace these paths we need to send another message
			current.times_visited += 1#blacklist ourselves for this one
			path.insert(0,current)
			current = current.pulse_pred[pulse_num]

		return path


"""
CENTRALIZED STUFF
"""


'''
use max flow on a subgraph of HG with all vertices at distance less than max_dist_scale * dist(s,t)

dist_measure ('path' or 't'): measure as closest distance to any vertex on the shortest path ('path') or as closest to t ('t')

This could be done decentralized according to the processes laid out in A. Segall's 1979 paper "decentralized maximum flow algorithms"

returns both the number of paths and the nodes used (assumed to be the entire subgraph)
'''
def hyper_VD_paths_local(HG:List[HyperNode],s:int,t:int,max_dist_scale=float('inf'),dist_measure='path',autoscale_increment=None):
	#first reduce HG down
	HGp = []
	nodes_removed = set()
	st_dist = hyper_dist(HG[s].coords,HG[t].coords)
	short_path = None
	if dist_measure == 'path':
		#find a shortest path from s to t
		short_path = HG[s].count_vd_paths_to_hyper(HG[t].coords,max_paths=1)[0]
	elif dist_measure == 't':
		short_path = [HG[t]]
	else:
		raise AttributeError('Unknown distance metric: {}'.format(dist_measure))

	for node in HG:
		if min([hyper_dist(node.coords,x.coords) for x in short_path]) <= max_dist_scale * st_dist:
			nodecpy = HyperNode(node.id,node.q,node.coords,node.idx,node.isom)
			nodecpy.neighbors = node.neighbors.copy()
			HGp.append(nodecpy)
		else:
			nodes_removed.add(node)

	#remove unnecessary edges
	for node in HGp:
		node.neighbors = node.neighbors - nodes_removed


	HGp_nx = convert_to_nx_graph(HGp)
	HGp_nx_transform = vertex_disjoint_transform(HGp_nx)

	if 'f{}'.format(s) in HGp_nx_transform.nodes:
		s_run = 'f{}'.format(s)
	else:
		s_run = s


	if autoscale_increment is not None:
		#keep increasing the max dist scale until the graph is connected
		ret = 0
		try:
			ret = nx.algorithms.flow.edmonds_karp(HGp_nx_transform,s_run,t).graph['flow_value'],len(HGp)
		except nx.NetworkXError:
			ret = hyper_VD_paths_local(HG,s,t,max_dist_scale=max_dist_scale+autoscale_increment,dist_measure=dist_measure,autoscale_increment=autoscale_increment)

		return ret
	else:
		ret = 0
		try:
			ret = nx.algorithms.flow.edmonds_karp(HGp_nx_transform,s_run,t).graph['flow_value'],len(HGp)
		except nx.NetworkXError:
			raise AttributeError('max dist scale of {} too small for s,t pair (s and t are cut) and autoscale is disabled.'.format(max_dist_scale))

		return ret