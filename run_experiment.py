"""
Framework for testing VD algorithms

this will always run some number of trials on some number of graphs, and output metrics in a csv format

#TODO ADD TESTING (UNIT TESTS)!!! (where?)
#TODO parallelize testing so it goes faster

#### Software engineering stuff

Function that takes in:
	- graph
	- algorithm to use
And outputs:
	- all these measures of how well it did

Script will: generate some graphs, run some algorithms on them, and output some numbers or plots
Example plot: GraphType of size 10, 100, 1000  (horizontal axis)
           vs Alg1 #paths_found, Alg2, ...
another for amount of compute (total? per node?)
(by the way: error bars)
another for number of nodes that need to participate?

In particular: ideally, a function that generates these graphs given number of desired nodes (and edges?)
	- "Erdos-Renyi"
	- Power-law or "social network" style
	- Grid sort of graph
	- Graphs we throw in for fun
"""
from TNUtil import *
from Hyperbolic import *
from stepper_sim import *
import argparse
import os.path

SYNCHRONIZED_PATH_ALGS = ['synch-v2']

DECENTRALIZED_PATH_ALGS = ['TN-v2',
						   'synch-v2',
						   'hyper-addr',
						   'hyper',
						   'hyper-neigh-addr',
						   'hyper-neigh',
						   'hyper-multi']

HYPER_EMBED_PATH_ALGS = ['hyper-addr',
						 'hyper',
						 'hyper-neigh-addr',
						 'hyper-neigh',
						 'hyper-multi']

CENTRALIZED_PATH_ALGS = ['only-opt','local-mf']

ALL_PATH_ALGS = DECENTRALIZED_PATH_ALGS + CENTRALIZED_PATH_ALGS

#inverse map from algorithms to the options they can use (that is, maps options to the algorithms that use them)
PATH_ALGS_OPTIONS_USED_INV = {'max_paths':['hyper','hyper-addr'],
							  'max_dist_scale':['hyper','hyper-addr','hyper-neigh','hyper-neigh-addr','hyper-multi'],
							  'stop_on_first_failure':['hyper','hyper-addr','hyper-neigh','hyper-neigh-addr','hyper-multi'],#TODO add this option for TN-v2?
							  #TODO add TTL
							  }

PATH_ALGS_NODE_TYPE = {x:TNNode for x in ALL_PATH_ALGS}#by default everything uses TNNodes
PATH_ALGS_NODE_TYPE.update({x:TNNode_Stepper for x in SYNCHRONIZED_PATH_ALGS})
PATH_ALGS_NODE_TYPE.update({x:HyperNode for x in HYPER_EMBED_PATH_ALGS})

'''
helper function for later things -- just pulling the code out of a later function here
'''
def validate_returned_paths(G:List[Union[TNNode,HyperNode,TNNode_Stepper]],paths:List[List[Union[int,TNNode,HyperNode,TNNode_Stepper]]],s:int,t:int,verbose=False):
	valid = True
	nodes_seen = set()
	for path in paths:
		prev_node = G[s]
		for node in path:
			# verify that paths exist
			if type(node) == int:
				node = G[node]
			if node not in prev_node.neighbors:
				valid = False
				if verbose:
					print('EDGE DOES NOT EXIST: ({},{})'.format(prev_node.id,node.id))
				else:
					break  # if we're printing everything we want to really pring EVERYTHING
			prev_node = node
			# verify that paths are disjoint
			if (node.id != s) and (node.id != t) and (node.id in nodes_seen):
				valid = False
				if verbose:
					print('REPEATED NODE ID: {}'.format(node.id))
				else:
					break  # if we're printing everything we want to really pring EVERYTHING
			else:
				nodes_seen.add(node.id)
		if prev_node.id != t:
			valid = False
			if verbose:
				print('DESTINATION NOT REACHED (FINAL NODE IS {}, t = {})'.format(prev_node,G[t]))

		if (not valid) and (not verbose):
			break
	return valid

'''
Run vd_path_alg on G from s to t

path alg possibilities:
	'TN-v2':	 		TNNode.count_vd_paths_to_v2(destination id, return paths = True)
	'synch-v2':			TNNode_Stepper.count_vd_paths_to_v2(destination id, return paths = True, TTL = infinity)
	'hyper-addr':		HyperNode.count_vd_paths_to_hyper_from_addr(dest address, npaths=inf, max distance scale = inf, stop on first fail = false)
	'hyper':			HyperNode.count_vd_paths_to_hyper(dest coordinates, npaths=inf, max distance scale = inf, stop on first fail = false)
	'hyper-neigh-addr':	HyperNode.count_vd_paths_to_hyper_neighborbl_from_addr(dest address, max distance scale = inf, stop on first fail = false)
	'hyper-neigh':		HyperNode.count_vd_paths_to_hyper_neighborbl(dest coordinates, max distance scale = inf, stop on first fail = false)
	'hyper-multi':		HyperNode.count_vd_paths_to_hyper_multivisit(dest coords, max dist scale = inf, stop on first fail = false)
	'local-mf':			hyper_VD_paths_local(hyperbolic graph, start, target, max distance scale = inf, distance measure = 'path', autoscale increment = none)
	'only-opt':			vertex_disjoint_paths(nx graph, start, target)
	(+more as they are added)


returns:
	success: bool: true iff all paths found were valid and all paths found were VD. if validate_paths=False, this is always true. if no paths are found, this is always true
	num_paths_found: how many paths did we find?
	num_nodes_used: how many nodes had to *send a message* -- nothing else counts towards this
	total_messages_sent: how many messages were sent?
	
does NOT return (i.e. these need to be calculated later):
	number of possible paths
	lower-bound network usage
	

this function will reset the graph after calculating the needed quantities
'''
def run_single_pair(G:List[Union[TNNode,HyperNode,TNNode_Stepper]],vd_path_alg:str,s:int,t:int,validate_paths=True,verbose=False,**kwargs):
	paths = []
	if vd_path_alg == 'only-opt':
		return True,0,0,0#basically, running just the optimal lets us skip all of this
	elif vd_path_alg == 'TN-v2':
		paths = G[s].count_vd_paths_to_v2(t,return_paths=True)
	elif vd_path_alg == 'synch-v2':
		if type(G) != List[TNNode_Stepper]:
			raise AttributeError("Synchronized algorithms can only be run on stepper-sim TNNodes")
		G[s].count_vd_paths_to_v2(t,return_paths=True,**kwargs)
		all_done_flag = False
		while not all_done_flag:
			all_done_flag = True
			for tn in G:
				tn.time += 1  # increment first so all of the nodes are incremented before we do this time step
			for tn in G:
				tn.increment_time()
			for tn in G:
				if len(tn.operations) > 0:  # has to be separate because high node ids can give operations to low ones
					all_done_flag = False

		G[s].cleanup()
		paths = G[s].paths
	elif vd_path_alg == 'hyper-addr':
		if type(G) != List[HyperNode]:
			raise AttributeError("Hyperbolic embedding algorithms can only be run on HyperNodes")
		paths = G[s].count_vd_paths_to_hyper_from_addr(G[t].saddr,**kwargs)
	elif vd_path_alg == 'hyper':
		if type(G[0]) != HyperNode:
			raise AttributeError("Hyperbolic embedding algorithms can only be run on HyperNodes")
		paths = G[s].count_vd_paths_to_hyper(G[t].coords,**kwargs)
	elif vd_path_alg == 'hyper-neigh-addr':
		if type(G[0]) != HyperNode:
			raise AttributeError("Hyperbolic embedding algorithms can only be run on HyperNodes")
		paths = G[s].count_vd_paths_to_hyper_neighborbl_from_addr(G[t].saddr,**kwargs)
	elif vd_path_alg == 'hyper-neigh':
		if type(G[0]) != HyperNode:
			raise AttributeError("Hyperbolic embedding algorithms can only be run on HyperNodes")
		paths = G[s].count_vd_paths_to_hyper_neighborbl(G[t].coords,**kwargs)
	elif vd_path_alg == 'hyper-multi':
		if type(G[0]) != HyperNode:
			raise AttributeError("Hyperbolic embedding algorithms can only be run on HyperNodes")
		paths = G[s].count_vd_paths_to_hyper_multivisit(G[t].coords,**kwargs)
	elif vd_path_alg == 'local-mf':
		if type(G[0]) != HyperNode:
			raise AttributeError("Hyperbolic embedding algorithms can only be run on HyperNodes")
		#we don't verify with this method since it doesn't return paths
		num_paths_found, num_nodes_used = hyper_VD_paths_local(G,s,t,**kwargs)
		return True,num_paths_found,num_nodes_used,0#message send calculation isn't done here since this is centralized
	#ADD NEW ALGORITHM PATH GENERATION SEMANTTICS HERE
	else:
		raise AttributeError("Unknown experiment algorithm: {}, see documentation for run_experiment::run_single_pair".format(vd_path_alg))

	valid = True
	if validate_paths:
		#make sure the paths are correct
		valid = validate_returned_paths(G,paths,s,t,verbose=verbose)


	#reset everything and tally up usage
	operations_this_pair = 0
	nodes_with_operations = 0
	for node in G:
		operations = node.reset_all_and_return_ops()
		if operations != 0:
			nodes_with_operations += 1
		operations_this_pair += operations

	return valid,len(paths),nodes_with_operations,operations_this_pair

'''
As with run_single_pair, but this runs many and calculates summary stats for the metrics

if sample_proportion is less than 1, we will randomly sample among the possible pairs of paths without replacement min(max_npairs,floor(num_possible_pairs*sample_proportion)) pairs to calculate

if progress interval is set at 0, no progress is printed, else progress wil be printed every <progress_interval> pairs

returns (mapping from metric string to tuple of mean and stdev)
	number of total paths found mean and standard deviation
	maximum possible paths found mean and stdev
	number of nodes used mean and standard deviation <approximate for centralized algorithms>
	optimal network usage (number of nodes) mean and standard deviation <only if "compare_to_opt" is true>
	messages sent per node (among those that send messages) mean and standard deviation <only for decentralized algorithms>
'''
def run_many_pairs(G:List[Union[TNNode,HyperNode,TNNode_Stepper]],vd_path_alg:str,sample_pair_proportion=1.,seed=None,max_npairs=float('inf'),progress_interval=0,compare_to_optimal=False,**kwargs):
	np.random.seed(seed)
	### SUMMARY STATS WILL BE CALCULATED ON THESE
	num_paths_found_agg = []
	num_opt_paths_agg = []
	num_nodes_used_agg = []
	optimal_usage_agg = []
	messages_sent_per_node_used_agg = []
	###

	#how many pairs could we do?
	num_edges = int(sum([len(u.neighbors) for u in G]) / 2)
	num_pairs_possible = int((len(G) * (len(G) - 1))/2) - num_edges #|V Ux V| - |E| = number of non-edge pairs in the graph (Ux is unordered cartesian product (size is half of (|V x V| - |V|), which will always be an even number) without self-pairs)

	#how many will we do?
	npairs = int(min(max_npairs,float(num_pairs_possible) * sample_pair_proportion))

	#what will the pairs be?
	pairs = pick_n_random_pairs(G,npairs)

	#do the pairs
	pairs_run = 0
	nxG = convert_to_nx_graph(G)
	for s,t in pairs:
		if (progress_interval != 0) and ((pairs_run % progress_interval) == 0):
			print('{} of {} pairs evaluated ({:.1f}%)'.format(pairs_run,npairs,100. * float(pairs_run) / float(npairs)))

		exact_total_paths = 0
		opt_num_nodes = 0
		if compare_to_optimal:
			exact_paths = vertex_disjoint_paths(nxG,s,t,retrace=True)
			if validate_paths:
				validate_returned_paths(G,exact_paths,s,t,verbose=True)
			opt_num_nodes = sum([len(path)-1 for path in exact_paths]) + 1#don't count t for all of these -- just count it once
			exact_total_paths = len(exact_paths)
		else:
			exact_total_paths = vertex_disjoint_paths(nxG,s,t,retrace=False)

		#get the numbers for this pair
		success,num_paths_found,num_nodes_used,num_messages_sent = run_single_pair(G,vd_path_alg,s,t,**kwargs)
		#TODO do something with the success variable or get rid of it

		#throw them in the aggregators
		num_paths_found_agg.append(num_paths_found)
		num_opt_paths_agg.append(exact_total_paths)
		num_nodes_used_agg.append(num_nodes_used)
		if compare_to_optimal:
			optimal_usage_agg.append(opt_num_nodes)
		if vd_path_alg not in CENTRALIZED_PATH_ALGS:
			if num_nodes_used != 0:
				messages_sent_per_node_used_agg.append(float(num_messages_sent)/float(num_nodes_used))
			else:
				messages_sent_per_node_used_agg.append(float('inf'))

		pairs_run += 1

	return npairs, num_paths_found_agg, num_opt_paths_agg, num_nodes_used_agg, optimal_usage_agg, messages_sent_per_node_used_agg

METRIC_ORDERING = ['num_paths_found_mean','num_paths_found_stdev',
					'num_opt_paths_mean','num_opt_paths_stdev',
					'num_nodes_used_mean','num_nodes_used_stdev',
					'optimal_usage_mean','optimal_usage_stdev',
					'messages_sent_per_node_used_mean','messages_sent_per_node_used_stdev']

GRAPH_FNS = {'con-er':generate_connected_ER_graph,  #connected erdos renyi
			 'pnas-sn':generate_connected_rand_graph_from_deg_dist,  #pnas social network
			 'con-er-var':generate_connected_variable_dense_ER_graph,  #connected erdos renyi -- variable density
			 'nbad-union':generate_nbad_unioning_fast,  #nbad example for path unioning vs best of returned paths (also for naive blacklisting)
			 'nbad-neighborbl':generate_nbad_neighborbl,  #nbad example for neighbor-blacklisting, assuming pathfinder is shortest-path (e.g. bfs)
			 'nbad-either':generate_nbad_either,  #nbad example for both unioning and neighbor-bl
			 #add map from graph type string onto the actual generating function here
			 }

GRAPH_TYPES = list(GRAPH_FNS.keys())

'''
Using the same graph generation algorithm, calculate all the numbers we want for many different sizes of graphs
'''
def run_many_pairs_on_many_random_graphs(graph_sizes,vd_path_alg:str,generator,show_progress=True,generator_args=(),generator_kwargs=None,node_args=(),node_kwargs=None,runner_kwargs=None):
	if node_kwargs is None:
		node_kwargs = {}
	if generator_kwargs is None:
		generator_kwargs = {}
	if runner_kwargs is None:
		runner_kwargs = {}

	results = []
	for ndone,size in enumerate(graph_sizes):
		if show_progress:
			print('Starting graph {} of {}'.format(ndone+1,len(graph_sizes)))
		#generate the graph
		generator_args = [size] + list(generator_args)#size must be dynamic, so it isn't already a part of the generator args
		G = generate_random_graph(generator,PATH_ALGS_NODE_TYPE[vd_path_alg],generator_args,generator_kwargs,node_args,node_kwargs)
		generator_args = generator_args[1:]#remove this size so that the next size calls the generator correctly
		rmp_ret = run_many_pairs(G,vd_path_alg,**runner_kwargs)

		results.append([len(G),*rmp_ret])
		if show_progress:
			print('Finished graph {} of {} ({:.1f}% done)'.format(ndone + 1,len(graph_sizes),100. * float(ndone + 1) / float(len(graph_sizes))))


	#deduce repetition number
	res_idx = 0
	while results[res_idx][0] == results[0][0]:
		res_idx += 1
	rep = res_idx

	#collect results by graph size
	res_idx = 0
	results_ret = []
	for size in graph_sizes[0:len(graph_sizes):rep]:
		#get all of the results that have this size
		this_size_results = results[res_idx:res_idx+rep]#these are the ones to be combined
		tmp_arr_res = [[] for _ in range(len(this_size_results[0])-2)]
		tmp_npairs = []
		for one_res in this_size_results:
			for i,arr_res in enumerate(one_res):
				if i >= 2:
					tmp_arr_res[i-2].extend(arr_res)
				elif i >= 1:
					tmp_npairs.append(arr_res)
		num_paths_found_agg,num_opt_paths_agg,num_nodes_used_agg,optimal_usage_agg,messages_sent_per_node_used_agg = tmp_arr_res
		# calculate summary stats
		num_paths_found_mean = np.mean(num_paths_found_agg)
		num_paths_found_stdev = np.std(num_paths_found_agg,ddof=1)
		num_opt_paths_mean = np.mean(num_opt_paths_agg)
		num_opt_paths_stdev = np.std(num_opt_paths_agg,ddof=1)
		num_nodes_used_mean = np.mean(num_nodes_used_agg)
		num_nodes_used_stdev = np.std(num_nodes_used_agg,ddof=1)
		tmp_res = {'num_paths_found_mean':num_paths_found_mean,
				   'num_paths_found_stdev':num_paths_found_stdev,
				   'num_opt_paths_mean':num_opt_paths_mean,
				   'num_opt_paths_stdev':num_opt_paths_stdev,
				   'num_nodes_used_mean':num_nodes_used_mean,
				   'num_nodes_used_stdev':num_nodes_used_stdev
				   }
		if compare_to_optimal:
			optimal_usage_mean = np.mean(optimal_usage_agg)
			optimal_usage_stdev = np.std(optimal_usage_agg,ddof=1)
			tmp_res.update({'optimal_usage_mean':optimal_usage_mean,
						'optimal_usage_stdev':optimal_usage_stdev})
		if vd_path_alg not in CENTRALIZED_PATH_ALGS:
			messages_sent_per_node_used_mean = np.mean(messages_sent_per_node_used_agg)
			messages_sent_per_node_used_stdev = np.std(messages_sent_per_node_used_agg,ddof=1)
			tmp_res.update({'messages_sent_per_node_used_mean':messages_sent_per_node_used_mean,
						'messages_sent_per_node_used_stdev':messages_sent_per_node_used_stdev})

		results_ret.append([size,np.mean(tmp_npairs),tmp_res])
		res_idx += rep

	return results_ret

"""
Given rough bounds, give me some graph sizes
"""
def generate_graph_sizes(num_sizes:int,min_size,max_size,num_repeat:int):
	if (num_sizes % num_repeat) != 0:
		print("WARNING: requested number of graphs ({}) is not compatible with number of size repetitions ({}).".format(num_sizes,num_repeat))
	graph_sizes_proto = np.linspace(min_size,max_size,num=num_sizes//num_repeat)
	#fill this with all the necessary repeats
	graph_sizes = []
	for gs in graph_sizes_proto:
		for _ in range(num_repeat):
			graph_sizes.append(int(np.round(gs)))

	return graph_sizes


#argument parsers yayayay

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Experiment runner for VD path algorithm comparison")
	parser.add_argument('-a','--path_algorithm',nargs=1,type=str,choices=ALL_PATH_ALGS,default=['only-opt'],help='Which algorithm will be used to calculate VD paths between s and t.')
	parser.add_argument('-G','--graph_type',nargs=1,type=str,choices=GRAPH_TYPES,default=[GRAPH_TYPES[0]],help='What type of graph (generator) should be used?')
	parser.add_argument('-0','--graph_arg_0',nargs=1,type=float,default=None,help="First graph generator argument (for ER graphs, this is the average degree, for directeds it's approximate reciprocity, for variable density, this is the scale)")
	parser.add_argument('-1','--graph_arg_1',nargs=1,type=str,default=None,help="Second graph generator argument (for variable density, this is the name of the function (see util::generate_connected_variable_dense_ER_graph))")
	parser.add_argument('-g','--num_graph_sizes',nargs=1,type=int,default=[10],help='How many different graph sizes will be used for testing')
	parser.add_argument('-l','--min_graph_size',nargs=1,type=int,default=[10],help='Minimum graph size to test')
	parser.add_argument('-b','--max_graph_size',nargs=1,type=int,default=[250],help='Maximum graph size to test')
	parser.add_argument('-r','--num_repeat',nargs=1,type=int,default=[1],help='How many times to repeat each graph size')
	parser.add_argument('-p','--sample_pair_proportion',nargs=1,type=float,default=[1.0],help='What proportion of the possible pairs should be tested?')
	parser.add_argument('-m','--max_npairs',nargs=1,type=float,default=[float('inf')],help='What is the maximum number of pairs to test for a given graph?')
	parser.add_argument('-S','--no_show_graph_progress',default=False,action='store_true',help='Should we not show progress at the graph-testing level? (by default, we will)')
	parser.add_argument('-i','--pair_progress_interval',nargs=1,type=int,default=[0],help="How often [after how many trials] should we show progress at the pair level? (default is zero; i.e. we dont' show progress at all)")
	parser.add_argument('-c','--optimal_usage',default=False,action='store_true',help='Should we calculate the optimal network usage?')
	parser.add_argument('-C','--validate_paths',default=False,action='store_true',help="Should we validate at every step that the paths returned are VD and correct?")
	parser.add_argument('-M','--max_paths',nargs=1,type=float,default=[float('inf')],help='(VD path parameter): what is the maximum number of paths to find between s and t?')
	parser.add_argument('-x','--max_dist_scale',nargs=1,type=float,default=[float('inf')],help='(Hyperbolic VD path parameter): what is the maximum distance from t as a multiple of the distance between s and t a node is allowed to be before it is not considered for routing?')
	parser.add_argument('-f','--stop_on_first_failure',default=False,action='store_true',help='(VD path parameter): should we stop looking for paths as soon as we fail to get a return from one neighbor? (can be a performance increase, but will also lose some paths)')
	parser.add_argument('-q','--hyper_embed_degree',nargs=1,type=int,default=[3],help="What degree should the hyperbolic embed addressing tree be? This must be one more than a power of two (i.e. 3,5,9,17...)")
	parser.add_argument('-o','--output',nargs=1,type=str,default=None,help="What file should we output to? (defaults to stdout)")
	parser.add_argument('-O','--overwrite_output_files',default=False,action='store_true',help="Should we overwrite preexisting output files? (default is to add a number at the end)")
	args = parser.parse_args()

	#reformat these in the way the program expects
	#graph size args begin
	num_graph_sizes = args.num_graph_sizes[0]
	min_graph_size = args.min_graph_size[0]
	max_graph_size = args.max_graph_size[0]
	num_repeat = args.num_repeat[0]
	#end
	#graph generator args begin
	generator = GRAPH_FNS[args.graph_type[0]]
	graph_generator_arguments = []
	for arg in [args.graph_arg_0
				,args.graph_arg_1,
				# add others here and in the parser (-1, -2, ...)
				]:
		if arg is not None:
			graph_generator_arguments.append(arg[0])
	#end
	#high level runner args begin
	sample_pair_proportion = args.sample_pair_proportion[0]
	max_npairs = args.max_npairs[0]
	show_graph_progress = not args.no_show_graph_progress
	#end
	#single pair args begin
	sp_progress = args.pair_progress_interval[0]
	compare_to_optimal = args.optimal_usage
	validate_paths = args.validate_paths
	output = args.output
	if output is not None:
		output = output[0]
	overwrite = args.overwrite_output_files
	#end
	#vd path alg args begin
	path_algorithm_str = args.path_algorithm[0]
	stop_on_first_failure = args.stop_on_first_failure
	path_algorithm_kwargs = {}
	for desc,val in zip(['max_paths','max_dist_scale','stop_on_first_failure'],[args.max_paths[0],args.max_dist_scale[0],stop_on_first_failure]):
		if path_algorithm_str in PATH_ALGS_OPTIONS_USED_INV[desc]:
			path_algorithm_kwargs.update({desc:val})
	# node args begin
	node_args = ()
	if path_algorithm_str in HYPER_EMBED_PATH_ALGS:
		q = args.hyper_embed_degree[0]
		node_args += (q,)
	# end

	#aggregate runner kwargs
	runner_kwargs = path_algorithm_kwargs
	runner_kwargs.update({'sample_pair_proportion':sample_pair_proportion,
						  'max_npairs':max_npairs,
						  'compare_to_optimal':compare_to_optimal,
						  'progress_interval':sp_progress,
						  })


	#prepare the things
	graph_sizes = generate_graph_sizes(num_graph_sizes,min_graph_size,max_graph_size,num_repeat)
	if (not overwrite) and (output is not None):
		#get the prefix and file extension (e.g. 'file' and '.csv' for output='file.csv')
		file_re_match = re.match(r'(.*)(\..+)',output)
		file_prefix = output
		file_type = ""
		if file_re_match is not None:
			file_prefix = file_re_match.group(1)
			file_type = file_re_match.group(2)

		#make sure the given output file name is okay (not already in use here)
		post_num = 0
		try_file = output
		while os.path.isfile(try_file):
			try_file = file_prefix + '_{}'.format(post_num) + file_type
			post_num += 1
		output = try_file

	#get the results
	results = run_many_pairs_on_many_random_graphs(graph_sizes,path_algorithm_str,generator,show_graph_progress,graph_generator_arguments,None,node_args,None,runner_kwargs)

	if output is None:
		#output to stdout
		print('\n')#get some space between the progress updater and the results
		#print the results
		print('graph_size,npairs',end='')
		#static metric ordering
		#print description first
		for metric in METRIC_ORDERING:
			if metric in results[0][2]:
				print(',{}'.format(metric,metric),end='')
		print()

		#then print values
		for graph_size,npairs,res_dict in results:
			print('{},{}'.format(graph_size,npairs),end='')
			for metric in METRIC_ORDERING:
				if metric in res_dict:
					print(',{}'.format(res_dict[metric]),end='')
			print()
	else:
		#output to outfile
		with open(output,'w') as out:
			# print the results
			out.write('graph_size,npairs')
			# static metric ordering
			# print description first
			for metric in METRIC_ORDERING:
				if metric in results[0][2]:
					out.write(',{}'.format(metric))
			out.write('\n')

			# then print values
			for graph_size,npairs,res_dict in results:
				out.write(str(graph_size))
				out.write(',{}'.format(npairs))
				for metric in METRIC_ORDERING:
					if metric in res_dict:
						out.write(',{}'.format(res_dict[metric]))
				out.write('\n')

		print('################################################## DONE ####################################################')

	#aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaand we're done. again. hopefully pycharm won't randomly decide to yeet half of this again