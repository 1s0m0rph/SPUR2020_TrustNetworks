#TODO remove this file when done with testing framework

from stepper_sim import *
from TN import *

#how many paths does this method find on average, as a percentage of the maximum (empirical)?

N = 100

np.random.seed(0)
tng = generate_rand_graph_from_deg_dist(N,approx_reciprocity=1,node_type=TNNode)

#check that all coords are unique
coords_seen = set()
for node in tng:
	tupcoords = tuple(node.coords)
	if tupcoords in coords_seen:
		print('COORDINATE DUPLICATED: {}, {}'.format(tupcoords,node))
	else:
		coords_seen.add(tupcoords)

# exit(0)#for now

#try now with TTL
init_TTL = float('inf')#we will consider paths of only this length or less

paths_found_sum = 0
paths_exact_sum = 0
num_network_used_sum = 0
node_operations_sum = 0
npairs = 0
pair_count = 0

for s in range(N):
	for t in range(s+1,N):
		if t not in tng[s].neighbors:
			npairs += 1#doing this first so we can get progress reports

# s = 3
# t = 8
# npairs = 1

for s in range(N):
	for t in range(s+1,N):
		if t not in tng[s].neighbors:
			if pair_count % 100 == 0:
				print('{} of {} ({:.3f}%)'.format(pair_count,npairs,100.*float(pair_count)/float(npairs)))
			exact_total_paths = vertex_disjoint_paths(convert_to_nx_graph(tng),s,t)

			paths = tng[s].count_vd_paths_to_v3(t,tng[t].coords,heuristic='dist')

			# print('{} of {} total paths found: '.format(len(paths),exact_total_paths))

			nodes_seen = set()

			for path in paths:
				prev_node = tng[s]
				for node in path:
					#verify that paths exist
					if node not in prev_node.neighbors:
						print('EDGE DOES NOT EXIST: ({},{})'.format(prev_node.id,node.id))
					prev_node = node
					#verify that paths are disjoint
					if (node.id != s) and (node.id != t) and (node.id in nodes_seen):
						print('REPEATED NODE ID: {}'.format(node.id))
					else:
						nodes_seen.add(node.id)
				# print(path)

			paths_found_sum += len(paths)
			paths_exact_sum += exact_total_paths
			pair_count += 1

			count_this_pair_network_used = 0#how many nodes had to do something this time around?
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
				tng[i].initial_TTL = float('inf')
				node_operations_sum += tng[i].operations_done
				if tng[i].operations_done != 0:
					count_this_pair_network_used += 1
				tng[i].operations_done = 0

			num_network_used_sum += count_this_pair_network_used

print('AVERAGE NUMBER (PCT) OF PATHS FOUND:\t\t\t\t\t{:.3f} ({:.3f}%)'.format(paths_found_sum/npairs, 100*(paths_found_sum/paths_exact_sum)))
print('AVERAGE NUMBER (PCT) OF NODES WITH OPERATIONS:\t\t\t{:.0f} ({:.3f}%)'.format(num_network_used_sum/npairs,100 * (num_network_used_sum / (npairs*len(tng)))))
print('AVERAGE NUMBER OF OPERATIONS PER NODE WITH OPERATIONS:\t{:.3f}'.format(node_operations_sum / num_network_used_sum))