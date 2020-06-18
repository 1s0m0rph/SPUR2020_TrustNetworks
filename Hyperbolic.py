from TN import *

'''
arguments are complex numbers
'''
def hyper_dist(a:complex,b:complex):
	return np.arccosh(1 + (2*(abs(a-b)**2))/((1 - abs(a)**2)*(1 - abs(b)**2)))


'''
convert a succinct address (in-person sharing) to coordinates by running the generators forward
'''
def addr_to_coords(q,addr):
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
		# for some reason doing it the way their pseudocode says to doesn't work because of this R^i thing
		# it only affects the zeroth generator (and therefore only the root)
		rot_isom_i = Isometry((np.array([rot_isom.rot,rot_isom.trans]) ** complex(i))[0],0j)
		generators.append(rot_isom_i.cross(trans_isom).cross(rot_isom_i.inv()))

	return generators

class HyperNode(TNNode):

	def __init__(self,node_id,q,coordinates=None,index=None,isometry=None,addr=None):
		super().__init__(node_id)
		self.coords = coordinates
		self.idx = index
		self.isom = isometry
		self.q = q

		self.d_coords = []  # list of available daughter coords (list of [(coords,index,isometry)] which can instantiate daughters)
		self.addr = addr#address
		# self.calculate_daughter_coords()

		self.neighbors = set()
		self.daughters = set()

		self.d_add_search_flag = False

	def __repr__(self):
		return "HN, {} (@{})".format(self.id,hex(self.addr)[2:])

	'''
	algorithm 2
	'''
	def init_as_root(self):
		self.coords = 0 + 0j
		self.idx = 0
		self.isom = Isometry(1 + 0j,0 + 0j)
		self.addr = 0
		self.calculate_daughter_coords()

	'''
	algorithm 3
	'''
	def calculate_daughter_coords(self):
		generators = define_generators(self.q)

		for i in range(0 if self.coords == 0j else 1,self.q):
			didx = (self.idx + i) % self.q
			disom = self.isom.cross(generators[didx])
			dcoord = disom.evaluate(0)
			w = int(np.log2(self.q+2))
			daddr = self.addr << w | (i + (1 if self.coords == 0j else 0))

			self.d_coords.append((dcoord,didx,disom,daddr))

	def reset_flags(self):
		super().reset_flags()
		self.d_add_search_flag = False

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
			# ask our neighbors if they can add d
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
			return  # already exists
		if self.coords is None:
			if n.coords is None:
				raise AttributeError("Tried to connect two nodes that weren't already in the network")
			# make n our parent
			self.coords, self.idx, self.isom, self.addr = n.add_daughter(self)
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
	def count_vd_paths_to_hyper_from_addr(self,dest_addr,npaths=float('inf'),max_dist_scale=float('inf')):
		dest_coords = addr_to_coords(self.q,dest_addr)
		return self.count_vd_paths_to_hyper(dest_coords,npaths=npaths,max_dist_scale=max_dist_scale)

	'''
	this should technically be private access
	'''
	def count_vd_paths_to_hyper(self,dest_coords,npaths=float('inf'),max_dist_scale=float('inf')):
		self.search_blacklist_flag = True
		st_dist = hyper_dist(self.coords,dest_coords)
		# start a search to dest from each neighbor
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

		for neighbor in self.neighbors:
			neighbor.reset_search()

		self.search_blacklist_flag = False

		return paths

	def gnh_interm(self,dest_coords,pred,st_dist,max_dist_scale):
		if self.coords == dest_coords:  # this has to happen before the visited check to implement path shadowing
			# blacklist nodes on this path
			self.pulse_pred.update({-1:pred})
			self.resetted_flag = False
			self.search_blacklist_flag = False  # t is never blacklisted, but the function after sets it to be so
			path = self.v2_vd_blacklist_zip(-1,[])
			return path

		if -1 in self.pulse_pred:
			return None  # we've already been visited in this search

		if self.search_blacklist_flag:
			return None  # we are already blacklisted (or have been visited in this search), so don't go this way

		# we've relayed this pulse now
		self.pulse_pred.update({-1:pred})
		self.resetted_flag = False

		# otherwise ask the *right* neighbor(s) if they know the muffin man
		neighbors_to_call = list(sorted(list(self.neighbors),key=lambda x:hyper_dist(x.coords,dest_coords)))
		for neighbor in neighbors_to_call:
			if hyper_dist(neighbor.coords,dest_coords) <= max_dist_scale * st_dist:
				self.operations_done += 1
				path_ret = neighbor.gnh_interm(dest_coords,self,st_dist,max_dist_scale)
				if path_ret is not None:
					return path_ret

		return None




def generate_rand_graph_from_deg_dist(num_nodes,q,dist=lambda:scipy.stats.truncnorm.rvs(0,float('inf'),loc=3,scale=3)):
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