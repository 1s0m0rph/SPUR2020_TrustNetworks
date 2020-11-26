"""
General abstraction on the idea of an embedding on a network

we should be able to apply this to TNNodes
"""
from abc import ABC, abstractmethod
from util import *
from TN import *



class Embedder(ABC):
	"""
	This is what we will use to annotate TNNodes
	"""

	EMBEDDINGS = ['hyperbolic','n-adic']

	def __init__(self,etype,atype,dtype):
		assert etype in self.EMBEDDINGS
		self.etype = etype#embedding type
		self.atype = atype#address type
		self.dtype = dtype#distance type

	def atype_check(self,obj):
		type_check(obj,self.atype)

	def dtype_check(self,obj):
		type_check(obj,self.dtype)

	def dist(self,adrx,adry):
		self.atype_check(adrx)
		self.atype_check(adry)
		d = self.__dist__(adrx,adry)
		self.dtype_check(d)
		return d

	@abstractmethod
	def __dist__(self,adrx,adry):
		pass

	@abstractmethod
	def address_graph(self,G:List[TNNode]):
		pass


class TreeEmbedder(Embedder,ABC):
	"""
	this is specifically for tree-like embeddings
	"""

	def __init__(self, etype, atype, dtype):
		super().__init__(etype,atype,dtype)

	def get_adr_children(self,adr):
		self.atype_check(adr)
		chls = self.__gac__(adr)
		for ch in chls:
			self.atype_check(ch)
		return chls.copy()

	@abstractmethod
	def __gac__(self,adr):
		pass

	def get_root(self):
		r = self.__grt__()
		self.atype_check(r)
		return r

	@abstractmethod
	def __grt__(self):
		pass

	def address_graph(self,G:List[TNNode],root_idx=None):
		if root_idx is None:
			root_idx = random.randint(0,len(G)-1)#random unless otherwise specified
		root_adr = self.get_root()
		#bfs on G to assign the addresses
		q = [G[root_idx]]
		G[root_idx].address = root_adr
		G[root_idx].adist = self.dist
		#we'll use the assignment of an address as a seen flag
		while len(q) > 0:
			current = q.pop(0)
			chl = self.__gac__(current.address)#don't be too picky about types with this one

			ai = 0
			for v in current.neighbors:
				if v.address is None:
					if ai >= len(chl):
						raise ValueError("The graph contains nodes of higher degree than the tree allows. Increase the tree branching parameter")
					v.address = chl[ai]
					v.adist = self.dist
					ai += 1
					q.append(v)


		for v in G:
			if (v.address is None) or (v.adist is None):
				raise ValueError("Graph is not connected! Embedding incomplete.")

		return G

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

class Hyperbolic_Embedder(TreeEmbedder):

	def __init__(self,q):
		super().__init__('hyperbolic',({complex,np.complex128},{int,np.int64,np.int32},Isometry),{float,np.float64})#TODO double check typing here
		self.q = q

	'''
	algorithm 3
	'''
	def calculate_daughter_coords(self,adr):
		coords,idx,isom = adr
		generators = define_generators(self.q)

		d_coords = []
		for i in range(self.q):
			didx = (idx+i)%self.q
			disom = isom.cross(generators[didx])
			dcoord = disom.evaluate(0)

			d_coords.append((dcoord,didx,disom))

		return d_coords

	def __dist__(self,adrx,adry):
		coordx,_,_ = adrx
		coordy,_,_ = adry
		return hyper_dist(coordx,coordy)

	def __gac__(self,adr):
		return self.calculate_daughter_coords(adr)

	def __grt__(self):
		rcoords = 0 + 0j
		ridx = 0
		risom = Isometry(1 + 0j,0 + 0j)
		return (rcoords,ridx,risom)



class NAdic(TreeEmbedder):

	def __init__(self,n):
		super().__init__("n-adic",n_adic,n_adic)
		if (type(n) != int) or (n < 2):
			raise TypeError("n-adic embedding requires n to be an integer greater than 1")
		self.n = n

	def __dist__(self,adrx,adry):
		return adrx.dist_to_DA(adry)

	def __gac__(self,adr):
		return adr.children()

	def __grt__(self):
		return n_adic(1,self.n,0,reduce=False,cache=False)


class n_adic:
	"""
	d = None is taken to mean d = -inf (denom = 0 or x = inf)
	"""


	'''
	cache can be any of 'all','chld','denom',(anything else results in no caching)
	'''
	def __init__(self, a: int, n: int, d, reduce=True, cache=False):
		assert ((d is None) or (type(d) == int))
		if reduce:
			while (a%n) == 0:
				a //= n
				d -= 1
		self.a = a
		self.n = n
		self.d = d
		self.cache = cache


		if d is None:
			self.chld = []
			self.truedenom = 0
		else:
			self.chld = None
			self.truedenom = None


	def __add__(self, other):
		if self.d > other.d:
			oscale = self.n**(self.d-other.d)
			sscale = 1
		elif other.d > self.d:
			oscale = 1
			sscale = self.n**(other.d-self.d)
		else:
			oscale = sscale = 1
		return n_adic(self.a*sscale+other.a*oscale, self.n, max(self.d, other.d), reduce=True,cache=self.cache)

	def __neg__(self):
		return n_adic(-self.a, self.n, self.d, reduce=False,cache=self.cache)

	def __sub__(self, other):
		return self+(-other)

	def __eq__(self, other):
		return (self.n == other.n) and (self.a == other.a) and (self.d == other.d)

	def __abs__(self):
		return n_adic(abs(self.a), self.n, self.d, reduce=False,cache=self.cache)

	def children(self):
		if self.chld is not None:
			return self.chld.copy()
		chdenom = self.d+1
		chld = []
		for chnum in range(self.n*(self.a-1)+1, self.n*self.a):
			chld.append(n_adic(chnum, self.n, chdenom, reduce=False,cache=self.cache))
		if self.cache:
			self.chld = chld
		return chld.copy()

	def is_ancestor_of(self, other):
		if self == other:
			return True
		if self.d >= other.d:
			return False
		#just need to check the main condition now
		scale = self.n**(other.d-self.d)
		rbound = scale*self.a
		if other.a >= rbound:
			return False
		lbound = rbound-scale
		return other.a >= lbound

	'''
	this version is ancestor-weighted ONLY (i.e. not descendant-weighted)
	'''

	def dist_to(self, other):
		raw_dist = abs(self-other)
		if self.is_ancestor_of(other):
			return raw_dist
		else:
			return raw_dist+n_adic(1, self.n, 0, reduce=False,cache=self.cache)

	'''
	this is the descendent-ancestor weighted metric
	'''
	def dist_to_DA(self, other):
		if self == other:  #both descendant and ancestor case
			return -n_adic(1, self.n, None, reduce=False,cache=self.cache)  #take this to mean (negative) infinity
		s_anc_o = self.is_ancestor_of(other)
		o_anc_s = other.is_ancestor_of(self)
		#we know that not both of these ^^ are true at this point
		if s_anc_o:
			return self.d-other.d  #-(e-d)
		elif o_anc_s:
			return other.d-self.d  #-(d-e)
		else:
			raw_dist = abs(self-other)
			return raw_dist

	def __repr__(self):
		if self.truedenom == 0:
			return ("-" if self.a < 0 else "")+"INF"
		tdenom = None
		if self.truedenom is not None:
			tdenom = self.truedenom
		elif self.cache:
			self.truedenom = self.n**self.d
		else:
			tdenom = self.n**self.d
		return "{}/{}".format(self.a, tdenom)


def lnat_inner(x, max_depth, depth=0):
	s = ''
	if depth > max_depth:
		return s

	s += '\t'*depth
	s += '[.'+str(x)
	s += '\n'

	#print all children
	for ch in x.children():
		s += lnat_inner(ch, max_depth, depth=depth+1)

	s += '\t'*depth
	s += ']\n'
	return s


def latex_n_adic_tree(n, depth):
	s = '\\Tree'
	x = n_adic(1, n, 0)
	s += lnat_inner(x, depth)
	print(s)


def get_embedder(etype: str):
	if etype == 'hyperbolic':
		return Hyperbolic_Embedder
	elif etype == 'n-adic':
		return n_adic