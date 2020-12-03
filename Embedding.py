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



class NAdic_Embedder(TreeEmbedder):

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
		return n_adic(1,self.n,0,n_adic_reduce=False)


def long_divide(a,b):
	q = a//b
	adiffr0 = a - q*b
	adiff0 = abs(adiffr0)
	adiffr1 = adiffr0 - b
	adiff1 = abs(adiffr1)
	if adiff0 < adiff1:
		return q,adiffr0
	else:
		return q+1,adiffr1

def ext_eucl_int(a:int,b:int,gcd_only=False):
	if a == 0:
		return b
	if b == 0:
		return a
	carda = (1,0)
	cardb = (0,1)

	q,r = long_divide(a,b)
	cardc = (carda[0] - (q*cardb[0]),carda[1] - (q*cardb[1]))
	carda = cardb
	cardb = cardc
	a = b
	b = r

	while r != 0:
		q, r = long_divide(a, b)
		cardc = (carda[0]-(q*cardb[0]), carda[1]-(q*cardb[1]))
		carda = cardb
		cardb = cardc
		a = b
		b = r

	if a < 0:
		a = -a
		carda = (-carda[0],-carda[1])

	if gcd_only:
		return a
	else:
		return a,carda

def gcd(a:int,b:int):
	return ext_eucl_int(a,b,gcd_only=True)

class Rational:
	"""
	a rational number a/b where a and b are coprime* integers (b = 0 is thought of as being infinite)
	
	* iff rat_reduce is set to true
	"""

	def __init__(self,a:int,b:int,rat_reduce=True):
		if rat_reduce:
			g = gcd(a,b)
			a //= g
			b //= g

		self.numerator = a
		self.b = b
		
	def pairwise_check(self,other):
		if type(other) in [int,np.int32,np.int64]:
			other = Rational(other,1,rat_reduce=False)
		elif type(other) in [float,np.float64]:
			other = continued_frac_approx_convergents(other)[-1]
		return other

	def __add__(self, other):
		other = self.pairwise_check(other)
		return Rational(self.numerator*other.b+other.numerator*self.b,self.b*other.b,rat_reduce=True)

	def __neg__(self):
		return Rational(-self.numerator,self.b)


	def __sub__(self, other):
		return self+(-other)

	def __mul__(self, other):
		other = self.pairwise_check(other)
		return Rational(self.numerator*other.numerator,self.b*other.b,rat_reduce=True)

	def __abs__(self):
		return Rational(abs(self.numerator),abs(self.b))

	def __repr__(self):
		if self.b == 0:
			return ("-" if self.numerator < 0 else "")+"INF"
		return "{}/{}".format(self.numerator,self.b)

	"""
	called with the unary "~" operator
	"""
	def __invert__(self):
		return Rational(self.b,self.numerator)

	def __truediv__(self, other):
		return self * (~other)

	def __floordiv__(self, other):
		return self / other

	def __le__(self, other):
		other = self.pairwise_check(other)
		diff = self - other
		return diff.numerator <= 0

	def __eq__(self, other):
		other = self.pairwise_check(other)
		return (self.numerator == other.numerator) and (self.b == other.b)

	def __lt__(self, other):
		return (self <= other) and (self != other)

	def __gt__(self, other):
		return not (self <= other)

	def __ge__(self, other):
		return not (self < other)


class n_adic(Rational):
	"""
	d = None is taken to mean d = -inf (denom = 0 or x = inf)
	"""


	def __init__(self,a:int,n: int,d,n_adic_reduce=True):
		assert ((d is None) or (type(d) == int))
		if n_adic_reduce:
			while ((a%n) == 0) and (d > 0):
				a //= n
				d -= 1
		self.numerator = a
		self.n = n
		self.exp = d

		if d is None:
			super().__init__(a,0,rat_reduce=False)
		else:
			super().__init__(a,n**d,rat_reduce=False)

	def __add__(self, other):
		other = self.pairwise_check(other)
		if self.exp is None:
			if (other.exp is None) and (((self.numerator < 0) and (other.numerator > 0)) or ((self.numerator > 0) and (other.numerator < 0))):
				raise ZeroDivisionError("INF + (-INF) is indeterminate")
			return self
		if other.exp is None:
			return other
		if self.exp > other.exp:
			oscale = self.n**(self.exp-other.exp)
			sscale = 1
		elif other.exp > self.exp:
			oscale = 1
			sscale = self.n**(other.exp-self.exp)
		else:
			oscale = sscale = 1
		return n_adic(self.numerator*sscale+other.numerator*oscale,self.n,max(self.exp,other.exp),n_adic_reduce=True)

	def __mul__(self, other):
		other = self.pairwise_check(other)
		#(a/n^d)*(b/n^e) = (ab)/(n^(d+e))
		return n_adic(self.numerator*other.numerator,self.n,self.exp+other.exp,n_adic_reduce=True)

	def __neg__(self):
		return n_adic(-self.numerator,self.n,self.exp,n_adic_reduce=False)

	def __eq__(self, other):
		return (self.n == other.n) and (self.numerator == other.numerator) and (self.exp == other.exp)

	def __abs__(self):
		return n_adic(abs(self.numerator),self.n,self.exp,n_adic_reduce=False)

	def pairwise_check(self,other):
		if type(other) in [int,np.int32,np.int64]:
			other = n_adic(int(other),self.n,0)
		elif type(other) in [float,np.float64]:
			other = continued_frac_nadic_approx(other,self.n)
		return other

	def children(self):
		chdenom = self.exp+1
		chld = []
		for chnum in range(self.n*(self.numerator-1)+1,self.n*self.numerator):
			chld.append(n_adic(chnum,self.n,chdenom,n_adic_reduce=False))
		return chld

	def is_ancestor_of(self, other):
		if self == other:
			return True
		if self.exp >= other.exp:
			return False
		#just need to check the main condition now
		scale = self.n**(other.exp-self.exp)
		rbound = scale*self.numerator
		if other.numerator >= rbound:
			return False
		lbound = rbound-scale
		return other.numerator >= lbound

	'''
	this version is ancestor-weighted ONLY (i.e. not descendant-weighted)
	'''
	def dist_to(self, other):
		raw_dist = abs(self-other)
		if self.is_ancestor_of(other):
			return raw_dist
		else:
			return raw_dist+n_adic(1,self.n,0,n_adic_reduce=False)

	'''
	this is the descendent-ancestor weighted metric
	'''
	def dist_to_DA(self, other):#FIXME what all about this fails to match the progress report? (which is correct)
		if self == other:  #both descendant and ancestor case
			return -n_adic(1,self.n,None,n_adic_reduce=False)  #take this to mean (negative) infinity
		s_anc_o = self.is_ancestor_of(other)
		o_anc_s = other.is_ancestor_of(self)
		#we know that not both of these ^^ are true at this point
		if s_anc_o:
			return n_adic(self.exp-other.exp,self.n,0,n_adic_reduce=False)  #-(e-d)#FIXME this part of the distance metric is wrong
		elif o_anc_s:
			return n_adic(other.exp-self.exp,self.n,0,n_adic_reduce=False)  #-(d-e)
		else:
			raw_dist = abs(self-other)
			return raw_dist


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

def eval_convergents(cidxs:List[int]):
	res = Rational(cidxs[-1],1)
	for cidx in cidxs[-2::-1]:
		res = ~res
		res += cidx
	return res

def continued_frac_convergents(r_inp:Rational) -> List[Rational]:
	#TODO is there a faster way than just brute forcing the actual convergents?
	r = abs(r_inp)
	i = r.numerator//r.b
	cidxs = [i]
	convs = [Rational(i,1)]
	rem = r - i
	while rem.numerator > 1:
		i = rem.b//rem.numerator
		rem = Rational(rem.b%rem.numerator,rem.numerator)
		cidxs.append(i)
		conv = eval_convergents(cidxs)
		convs.append(conv)
	convs.append(r)
	return convs

def continued_frac_approx_convergents(x:Union[float,np.float64],w=100) -> List[Rational]:
	if not np.isfinite(x):
		return [Rational(int(np.sign(x)),0)]
	#first generate a totally brain-dead guess (i.e. <integer part of x> + <rational part of x>*2^w / 2^w
	i = int(x)
	ratxnum = int((x-i)*(2**w))
	if ratxnum == 0:
		return [Rational(i,1)]
	rat = Rational(ratxnum,1<<w) + i
	convs = continued_frac_convergents(rat)
	return convs

'''
w is the truncation width for the rational part of x
'''
def continued_frac_nadic_approx(x:Union[float,np.float64],n:int,w=100) -> n_adic:
	convs = continued_frac_approx_convergents(x,w)
	if convs[0].b == 1:
		if len(convs) != 1:
			convs = convs[1:]#drop things like /1 and what not
	if convs[0].b == 0:
		return n_adic(convs[0].numerator,n,0)
	if len(convs) == 1:
		return n_adic(convs[0].numerator,n,int(np.round(np.log(convs[0].b)/np.log(n))))


	#pick the closest one to some a/n^d with d < max_n_exp
	#this amounts to picking the one whose denom is closest to a power of n
	#which amounts to picking the one which has log_n(denom) closest to an integer (which will be d)
	lscale = 1/np.log(n)
	res_opt = np.log(convs[0].b)*lscale
	resd = int(np.round(res_opt))
	resn = convs[0].numerator
	res_opt = abs(res_opt - resd)
	for conv in convs[1:]:
		opt = np.log(conv.b)*lscale
		vd = int(np.round(opt))
		vn = conv.numerator
		opt = abs(opt - vd)
		if opt <= res_opt:
			resd = vd
			resn = vn
			res_opt = opt

	return n_adic(resn,n,resd)


def latex_n_adic_tree(n, depth):
	s = '\\Tree'
	x = n_adic(1, n, 0)
	s += lnat_inner(x, depth)
	print(s)


def get_embedder(etype: str):
	if etype == 'hyperbolic':
		return Hyperbolic_Embedder
	elif etype == 'n-adic':
		return NAdic_Embedder