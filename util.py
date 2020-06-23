import numpy as np
import networkx as nx
import rsa
import scipy.stats

'''
'add one' (LE) to the set (thinking of it is a bitvector/binary set)

universe is a list to maintain ordering
'''
def set_increment(universe:list,current:set):
	for elt in universe:
		if elt not in current:
			current.add(elt)
			return current
		else:
			current.remove(elt)

	return current