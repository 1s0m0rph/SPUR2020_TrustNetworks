#TODO remove this file when done with testing framework

from Hyperbolic import *

N = 100
q = int((2**1) + 1)

np.random.seed(0)
tng = generate_rand_graph_from_deg_dist_hyper(N,q)
max_dist_scale = 1
dist_measure = 'path'
autoscale_increment = 0.1

#verify uniqueness of coordinates and addresses
coords_seen = set()
addrs_seen = set()
for node in tng:
	if node.coords in coords_seen:
		print('COORDINATE DUPLICATED: {}, {}'.format(node.coords,node))
	elif node.saddr in addrs_seen:
		print('ADDRESS DUPLICATED: {}, {}'.format(node.saddr,node))
	else:
		coords_seen.add(node.coords)
		addrs_seen.add(node.saddr)

# #test succinct functions
# for node in tng:
# 	res = addr_to_coords(q,node.saddr,node.ADDRESS_LEVEL_BITS)
# 	if not np.isclose(res,node.coords):
# 		print('SUCCINCT ADDRESS FAILURE: {}, {}, {}'.format(node.coords,res,node.saddr))
#
# for node in tng:
# 	print(node)
#
# exit(0)

paths_found_sum = 0
paths_exact_sum = 0
num_network_used_sum = 0
node_operations_sum = 0
npairs = 0
pair_count = 0

for s in range(N):
	for t in range(s+1,N):
		if tng[t] not in tng[s].neighbors:
			npairs += 1#doing this first so we can get progress reports

# s = 0
# t = 4
# npairs = 1

for s in range(N):
	for t in range(s+1,N):
		if tng[t] not in tng[s].neighbors:
			if pair_count % 100 == 0:
				print('{} of {} ({:.3f}%)'.format(pair_count,npairs,100.*float(pair_count)/float(npairs)))
			exact_total_paths = vertex_disjoint_paths(convert_to_nx_graph(tng),s,t)

			num_paths_found,num_nodes_used = hyper_VD_paths_local(tng,s,t,max_dist_scale=max_dist_scale,dist_measure=dist_measure,autoscale_increment=autoscale_increment)
			# print('{} of {} total paths found: '.format(len(paths),exact_total_paths))


			paths_found_sum += num_paths_found
			paths_exact_sum += exact_total_paths
			pair_count += 1

			num_network_used_sum += num_nodes_used

print()
print('AVERAGE NUMBER (PCT) OF PATHS FOUND:\t\t\t\t\t{:.3f} ({:.3f}%)'.format(paths_found_sum/npairs, 100*(paths_found_sum/paths_exact_sum)))
print('AVERAGE NUMBER (PCT) OF NODES WITH OPERATIONS:\t\t\t{:.0f} ({:.3f}%)'.format(num_network_used_sum/npairs,100 * (num_network_used_sum / (npairs*len(tng)))))