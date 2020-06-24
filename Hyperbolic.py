from TN import *

'''
arguments are complex numbers
'''
def hyper_dist(a:complex,b:complex):
	return np.arccosh(1 + (2*(abs(a-b)**2))/((1 - abs(a)**2)*(1 - abs(b)**2)))

'''
modified to work with hyperbolic-embed graphs
'''
def generate_rand_graph_from_deg_dist_hyper(num_nodes,q,dist=lambda:scipy.stats.truncnorm.rvs(0,float('inf'),loc=3,scale=3)):
	#use the social network rand gen algorithm to make a new network
	G = [HyperNode(i,q) for i in range(num_nodes)]
	#root is 0
	G[0].init_as_root()

	#assign a degree to each node
	degrees = [max(1,int(dist())) for i in range(num_nodes)]#make sure everything is connected to at least one
	connections_left_to_make = {i for i in range(num_nodes)}#these nodes all have connections left to make
	total_connections_to_make = sum(degrees)
	nodes_in_network = {0}#only the root is in the network to start with

	#randomly connect the nodes according to their degrees
	connections_made = 0
	while (connections_made < total_connections_to_make) and (len(connections_left_to_make) > 1):
		#pick a random starter node (among the nodes in the network already
		i = np.random.choice(list(connections_left_to_make.intersection(nodes_in_network)))
		#pick a random other node to assign to
		assign_to_idx = np.random.choice(list(connections_left_to_make - {i} - G[i].neighbors))#no self-edges, no multiedges
		assign_to = G[assign_to_idx]
		nodes_in_network.add(assign_to_idx)#this is why we use a set for nodes in network
		G[i].add_public_key_in_person(assign_to)
		connections_made += 2
		degrees[i] -= 1
		if degrees[i] == 0:
			connections_left_to_make.remove(i)
		degrees[assign_to_idx] -= 1
		if degrees[assign_to_idx] == 0:
			connections_left_to_make.remove(assign_to_idx)

		#reciprocity is guaranteed by the add_neighbor method

	return G

'''
undoes the compression from str_compress (note: outputs a string of hex)
'''
def str_decompress(x:str):
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


'''
convert a succinct address (in-person sharing) to coordinates by running the generators forward
'''
def addr_to_coords_v0(q,addr):
	#TODO is there a more succinct way to do this? With 1000 nodes and q=30 we're getting >64bit addresses
	dummy = HyperNode(-1,q)
	dummy.init_as_root()
	addr_str = bin(addr)[2:]#should be unsigned so this shouldn't matter
	#pad with zeroes until the length is correct
	w = int(np.log2(q + 2))#q should be 2 less than a power of 2
	while len(addr_str) % w != 0:
		addr_str = '0' + addr_str

	gen_idx = int(addr_str[:w],2) - 1#generator index is funky for the first one so that root has a unique address
	if gen_idx == -1:
		return dummy.coords
	dummy = HyperNode(-1,q,*dummy.d_coords[gen_idx])
	dummy.calculate_daughter_coords()
	addr_str = addr_str[w:]
	while len(addr_str) > 0:
		gen_idx = int(addr_str[:w],2) - 1
		if gen_idx == -1:
			return dummy.coords
		dummy = HyperNode(-1,q,*dummy.d_coords[gen_idx])
		dummy.calculate_daughter_coords()
		addr_str = addr_str[w:]

	return dummy.coords


'''
convert a succinct address (in-person sharing) to coordinates by running the generators forward

version 1
'''
def addr_to_coords(q,addr,n_levelbits):
	dummy = HyperNode(-1,q)
	dummy.ADDRESS_LEVEL_BITS = n_levelbits
	dummy.init_as_root()
	addr = dummy.succ_addr_to_int(addr)
	addr_str = bin(addr)[2:]#should be unsigned so this shouldn't matter
	#pad with zeroes until the length is correct
	w = int(np.log2(q-1))#q should be 1 more than a power of 2
	while len(addr_str) < np.dtype(addr).itemsize*8:
		addr_str = '0' + addr_str

	#get the level
	level = int(addr_str[:n_levelbits],2)
	addr_str = addr_str[n_levelbits:]

	gen_idx = 0
	addr_str = (''.join(reversed(addr_str)))[:level*w]#have to go backwards because of how things were added
	for _ in range(level):
		gen_idx = int(addr_str[-w:],2)
		dummy = HyperNode(-1,q,*dummy.d_coords[gen_idx])
		dummy.calculate_daughter_coords()
		addr_str = addr_str[:-w]

	return dummy.coords


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

	ADDRESS_DTYPE = np.int64
	ADDRESS_LEVEL_BITS = 8#how many bits in the address are devoted to the level in the tree?

	def __init__(self,node_id,q,coordinates=None,index=None,isometry=None,saddr=None):
		super().__init__(node_id)
		self.coords = coordinates
		self.idx = index
		self.isom = isometry
		self.q = q

		self.d_coords = []#list of available daughter coords (list of [(coords,index,isometry)] which can instantiate daughters)
		self.saddr = saddr#address (k bits of level, width-k bits of address)
		#only store the succinct address since the coordinates are a better exact address anyway
		#self.calculate_daughter_coords()

		self.neighbors = set()
		self.daughters = set()

		self.d_add_search_flag = False

		self.pulse_blacklisted = set()#which pulses am I blacklisted for?
		self.max_neighbor_called = -1

	def __repr__(self):
		return "HN, {} (@{})".format(self.id,self.saddr)


	'''
	convert an int address to a succinct address
	'''
	def int_addr_to_succ(self,iaddr):
		n_nybble = np.dtype(self.ADDRESS_DTYPE).itemsize * 2
		uncompressed = hex(iaddr)[2:]
		while len(uncompressed) < n_nybble:
			uncompressed = '0' + uncompressed
		return str_compress(uncompressed)

	'''
	convert a succinct address to an int address
	'''
	def succ_addr_to_int(self,saddr):
		return self.ADDRESS_DTYPE(int(str_decompress(saddr),16))


	'''
	algorithm 2
	'''
	def init_as_root(self):
		self.coords = 0 + 0j
		self.idx = 0
		self.isom = Isometry(1 + 0j,0 + 0j)
		self.saddr = '0'#reserved root address
		self.calculate_daughter_coords()

	'''
	algorithm 3
	'''
	def calculate_daughter_coords(self):
		generators = define_generators(self.q)

		paddr = self.succ_addr_to_int(self.saddr)

		w = self.ADDRESS_DTYPE(int(np.log2(self.q - 1)))
		level_mask = (self.ADDRESS_DTYPE(1) << (np.dtype(self.ADDRESS_DTYPE).itemsize * 8 - 1)) >> self.ADDRESS_LEVEL_BITS - 1
		d_level = (((paddr & level_mask) >> (np.dtype(self.ADDRESS_DTYPE).itemsize * 8 - self.ADDRESS_LEVEL_BITS)) + 1)
		bound_check_mask = ~((self.ADDRESS_DTYPE(1) << (np.dtype(self.ADDRESS_DTYPE).itemsize * 8 - 1)) >> (np.dtype(self.ADDRESS_DTYPE).itemsize * 8 - self.ADDRESS_LEVEL_BITS))
		if d_level >= bound_check_mask:
			raise OverflowError('{} level bits insufficient to properly index tree'.format(self.ADDRESS_LEVEL_BITS))
		d_level <<= (np.dtype(self.ADDRESS_DTYPE).itemsize * 8 - self.ADDRESS_LEVEL_BITS)



		for i in range(1,self.q):#ignore the zeroth daughter of the root to enforce regularity of addressing
			didx = (self.idx + i) % self.q
			disom = self.isom.cross(generators[didx])
			dcoord = disom.evaluate(0)
			daddr = ((paddr if self.coords != 0j else self.ADDRESS_DTYPE(0)) & ~level_mask) << w | self.ADDRESS_DTYPE(i - 1)
			daddr |= d_level
			daddr = self.int_addr_to_succ(daddr)

			self.d_coords.append((dcoord,didx,disom,daddr))

	def reset_flags(self):
		super().reset_flags()
		self.d_add_search_flag = False
		self.pulse_blacklisted = set()
		self.max_neighbor_called = -1

	'''
	add this node as a daughter and give it the info it needs (coords, index, isometry)
	'''
	def add_daughter(self,d):
		if self.d_add_search_flag:
			return None
		self.d_add_search_flag = True
		self.resetted_flag = False
		if len(self.d_coords) > 0:
			self.daughters.add(d)
			return self.d_coords.pop(0)
		else:
			#ask our neighbors if they can add d
			for neighbor in self.neighbors:
				info = neighbor.add_daughter(d)
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
			self.coords, self.idx, self.isom, self.saddr = n.add_daughter(self)
			n.reset_search()
			self.calculate_daughter_coords()

		self.neighbors.add(n)
		n.add_neighbor(self)#enforce reciprocity

	def add_public_key_in_person(self,t_node):
		self.add_neighbor(t_node)


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
	def count_vd_paths_to_hyper_from_addr(self,dest_addr,npaths=float('inf'),max_dist_scale=float('inf'),stop_on_first_failure=False):
		dest_coords = addr_to_coords(self.q,dest_addr,self.ADDRESS_LEVEL_BITS)
		return self.count_vd_paths_to_hyper(dest_coords,npaths=npaths,max_dist_scale=max_dist_scale,stop_on_first_failure=stop_on_first_failure)

	'''
	this should technically be private access
	'''
	def count_vd_paths_to_hyper(self,dest_coords,npaths=float('inf'),max_dist_scale=float('inf'),stop_on_first_failure=False):
		self.search_blacklist_flag = True
		st_dist = hyper_dist(self.coords,dest_coords)
		#start a search to dest from each neighbor
		neighbors_to_call = list(sorted(list(self.neighbors),key=lambda x: hyper_dist(x.coords,dest_coords)))
		paths = []
		for neighbor in neighbors_to_call:
			if hyper_dist(neighbor.coords,dest_coords) <= max_dist_scale * st_dist:
				self.operations_done += 1
				path_ret = neighbor.gnh_interm(dest_coords,self,st_dist,max_dist_scale)
				if path_ret is not None:
					paths.append(path_ret)
					if len(paths) >= npaths:
						break
				elif stop_on_first_failure:
					break

		for neighbor in self.neighbors:
			neighbor.reset_search()

		self.search_blacklist_flag = False

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
		neighbors_to_call = list(sorted(list(self.neighbors),key=lambda x:hyper_dist(x.coords,dest_coords)))
		for neighbor in neighbors_to_call:
			if hyper_dist(neighbor.coords,dest_coords) <= max_dist_scale * st_dist:
				self.operations_done += 1
				path_ret = neighbor.gnh_interm(dest_coords,self,st_dist,max_dist_scale)
				if path_ret is not None:
					return path_ret

		return None
	
	
	"""
	partial/multiple blacklisting (doesn't seem to really work at all)
	"""

	'''
	initialize the search

	npaths variable tells us how many paths to find (we'll stop when we find this many or when we have found the max).
	max distance scale tells us how far nodes are allowed to be from t, as a linear function of the distance between s and t
		specifically, nodes that are further than max_dist_scale * (dist from s to t) are excluded
	'''
	def count_vd_paths_to_hyper_multibl_from_addr(self,dest_addr,npaths=float('inf'),max_dist_scale=float('inf')):
		dest_coords = addr_to_coords(self.q,dest_addr,self.ADDRESS_LEVEL_BITS)
		return self.count_vd_paths_to_hyper_multibl(dest_coords,npaths=npaths,max_dist_scale=max_dist_scale)


	'''
	this should technically be private access
	'''
	def count_vd_paths_to_hyper_multibl(self,dest_coords,npaths=float('inf'),max_dist_scale=float('inf')):
		st_dist = hyper_dist(self.coords,dest_coords)
		#start a search to dest from ourselves
		candidate_paths = []
		pulse_num = 0
		self.operations_done += 1
		path_ret = self.gnh_multibl_interm(dest_coords,None,pulse_num,st_dist,max_dist_scale)
		pulse_num = 1#don't send the 0-pulse twice
		while path_ret is not None:
			candidate_paths.append(path_ret)
			path_ret = self.gnh_multibl_interm(dest_coords,None,pulse_num,st_dist,max_dist_scale)
			self.operations_done +=  1
			pulse_num +=  1

		for neighbor in self.neighbors:
			neighbor.reset_search()

		#now pick the maximum simultaneous paths from the candidates

		compat_table = [[True for _ in range(len(candidate_paths))] for _ in range(len(candidate_paths))]#which paths are compatible with which others?
		#build the compatibility table
		for i,ipath in enumerate(candidate_paths):
			for j,jpath in enumerate(candidate_paths):
				#check these paths for compatibility
				seen_nodes_i = set()
				seen_nodes_j = set()
				for k in range(min(len(ipath)-1,len(jpath)-1)):#don't count the last one no matter what (it's t)
					if (ipath[k] == jpath[k]) or (ipath[k] in seen_nodes_j) or (jpath[k] in seen_nodes_i):
						compat_table[i][j] = False
						break
					else:
						seen_nodes_i.add(ipath[k])
						seen_nodes_j.add(jpath[k])

				#verify the longer one against the existing set
				lpath = ipath if len(ipath) > len(jpath) else jpath
				for k in range(min(len(ipath)-1,len(jpath)-1),len(lpath)-1):
					if (lpath[k] in seen_nodes_i) or (lpath[k] in seen_nodes_j):#either we have an invalid path or these are just incompatible
						compat_table[i][j] = False
						break

		#evaluate the compatibility table
		path_idxs = list(range(len(candidate_paths)))
		current_use_set = set_increment(path_idxs,set())#which paths (indexes) are we currently using?
		max_use_set = set()#which paths (indexes) form the best set so far?
		while len(current_use_set) > 0:
			if len(current_use_set) > len(max_use_set):#speedup: don't bother if it wouldn't be better anyway
				#check if this set is mutually compatible
				is_mutual = True
				for i in current_use_set:
					for j in current_use_set:
						if j < i:
							if not compat_table[i][j]:
								is_mutual = False#there is some pair that doesn't work
								break
					if not is_mutual:
						break

				if is_mutual:
					max_use_set = current_use_set.copy()
			current_use_set = set_increment(path_idxs,current_use_set)


		return [candidate_paths[i] for i in max_use_set]


	def gnh_multibl_interm(self,dest_coords,pred,pulse_num,st_dist,max_dist_scale):
		if self.coords == dest_coords:#this has to happen before the visited check to implement path shadowing
			#blacklist nodes on this path
			self.pulse_pred.update({pulse_num:pred})
			self.resetted_flag = False
			path = self.multi_blacklist_zip(pulse_num,[])
			return path

		if pulse_num in self.pulse_pred:
			return None#we've already been visited in this search

		if pulse_num in self.pulse_blacklisted:
			return None#we are already blacklisted (or have been visited in this search), so don't go this way

		#we've relayed this pulse now
		self.pulse_pred.update({pulse_num:pred})
		self.resetted_flag = False

		#otherwise ask the *right* neighbor(s) if they know the muffin man
		neighbors_to_call = list(sorted(list(self.neighbors),key=lambda x:(len(x.pulse_blacklisted),hyper_dist(x.coords,dest_coords))))#prioritize unblacklisted neighbors
		start_idx = 0 if (len(self.pulse_blacklisted) == 0) else (self.max_neighbor_called + 1)
		for nidx in range(start_idx,len(neighbors_to_call)):
			neighbor = neighbors_to_call[nidx]
			if hyper_dist(neighbor.coords,dest_coords) <= max_dist_scale * st_dist:
				self.operations_done += 1
				path_ret = neighbor.gnh_multibl_interm(dest_coords,self,pulse_num,st_dist,max_dist_scale)
				self.max_neighbor_called = nidx
				if path_ret is not None:
					return path_ret

		return None

	'''
	I am on the path to dest from s, so I need to be blacklisted

	this method also reconstructs the path
	'''
	def multi_blacklist_zip(self,pulse_num,path: list):
		if (pulse_num not in self.pulse_pred) or (self.pulse_pred[pulse_num] is None):
			return path #drop the first one (s)

		self.pulse_blacklisted.add(pulse_num)
		return self.pulse_pred[pulse_num].multi_blacklist_zip(pulse_num,[self] + path)


"""
CENTRALIZED STUFF
"""


'''
use max flow on a subgraph of HG with all vertices at distance less than max_dist_scale * dist(s,t)

dist_measure ('path' or 't'): measure as closest distance to any vertex on the shortest path ('path') or as closest to t ('t')

This could be done decentralized according to the processes laid out in A. Segall's 1979 paper "decentralized maximum flow algorithms"

returns both the number of paths and the nodes used (assumed to be the entire subgraph)
'''
def hyper_VD_paths_local(HG:list,s:int,t:int,max_dist_scale=float('inf'),dist_measure='path',autoscale_increment=None):
	#first reduce HG down
	HGp = []
	nodes_removed = set()
	st_dist = hyper_dist(HG[s].coords,HG[t].coords)
	short_path = None
	if dist_measure == 'path':
		#find a shortest path from s to t
		short_path = HG[s].count_vd_paths_to_hyper(HG[t].coords,npaths=1)[0]
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