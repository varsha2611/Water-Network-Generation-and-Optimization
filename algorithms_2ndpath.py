# -*- coding: utf-8 -*-
'''
Multiscale Entropic Network Generator (MUSKETEER)

Copyright (c) 2011-2012 by Alexander Gutfraind and Ilya Safro. 
All rights reserved.

Use and redistribution of this file is governed by the license terms in
the LICENSE file found in the project's top-level directory.

Core Algorithms 

wishlist
    the parameter 'level' should be eliminated for aesthetics of the recursion
        we could have helper method level_data, deeper_data = retrieve_level_data(params)
        retrieve_level_data should not only pull out the data but also do a sanity check
    support high-level control of editing: "fine", "medium, "coarse"
        (using recursion)
        

'''

import os
import time
import numpy as np
import numpy.random as npr
import random, sys
import networkx as nx
import pdb
import pickle
import gc
import simpletesters
import graphutils
import alternatives #module might be refered in params['algorithm']
import collections
try:
    import new_algs #for testing
except:
    pass
np.seterr(all='raise')

timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())
random_one = lambda arr: random.sample(arr, 1)[0]
max_int = np.iinfo(np.int32(10)).max

def chart_paths(G, source, target_node, new_edge_horizon, search_method, params):
#determine the path structure for paths from source to target_node
    estimate_of_paths = np.zeros(new_edge_horizon+1)
    num_mutual_nbs    = None #for nodes at distance d=2 only
    if search_method == 'bypass_path_only':  #find the length of just the bypass path = the path of second-shortest length
        distance = set()
        fringe = G.neighbors(source)
        fringe.remove(target_node)  #block the direct path
        target_found = False
        for d in range(1, new_edge_horizon):
            new_fringe = []
            for v in fringe:
                if v == target_node:
                    target_found = True
                    break
                if v not in distance:
                    distance.add(v)
                    new_fringe += G.neighbors(v)
            fringe = new_fringe
            if target_found:
                break
        if (not target_found) and (target_node in new_fringe):
            target_found = True
            d += 1
        if target_found:
            estimate_of_paths[d] = 1.0
    elif search_method == 'flood': #plot all the routes to target_node
        fringe = G.neighbors(source)
        fringe.remove(target_node)  #block the direct path
        target_found = False
        for d in range(1, new_edge_horizon):
            new_fringe = []
            for v in fringe:
                if v == source:
                    continue
                if v == target_node:
                    estimate_of_paths[d] += 1
                    target_found = True
                    continue
                new_fringe += G.neighbors(v)
            fringe = new_fringe 
            #each node in new_fringe represents a route;  
            #only routes that stay in one place are not counted; 
            #backtracking is allowed, which is unfortunate
            #(to exclude we would need to explicitly remember each route)
        if (not target_found) and (target_node in new_fringe):
            target_found = True
            #correct but unnecessary: 
            #estimate_of_paths[d+1] = new_fringe.count(target_node) 
            estimate_of_paths[d+1] = 1.0 
        if target_found: 
            estimate_of_paths /= sum(estimate_of_paths) 
    elif search_method == 'particles':
        #we send 100 particles on random paths, making sure that no particle is making loops
        # and for each d, measure how many of the particles reach the target on path of distance d.
        # potentially, we could relax the horizon somewhat b/c the number of particles decreasing rapidly
        good_nbs = G.neighbors(source)
        good_nbs.remove(target_node)  #block the direct path
        particles = [[random.choice(good_nbs)] for i in range(100)] #start location of the particles
        #todo: it would be more efficient to store in each particle its current location and a set of all previous locations (faster lookup)
        target_found = False
        for d in range(2, new_edge_horizon+1):
            for v in particles:
                last_loc = v[-1]
                next_loc = random.choice(G.neighbors(last_loc))
                if next_loc == source or next_loc in v:
                    particles.remove(v)
                elif next_loc == target_node:
                    estimate_of_paths[d] += 1
                    target_found = True
                    particles.remove(v)
                else:
                    v.append(next_loc)
        if target_found: 
            estimate_of_paths /= sum(estimate_of_paths) 
    else:
        raise ValueError('Unknown search method')
    
    if estimate_of_paths[2] > 0:
        num_mutual_nbs = len(set(G.neighbors(source)).intersection(G.neighbors(target_node)))
    return estimate_of_paths, target_found, num_mutual_nbs
        

def check_and_fix_connectivity(G, params):
#wishlist: this code could also consider tpl_data and insert so as not to affect clustering etc
    new_edges = set() 
    ccs = nx.connected_components(G) 
    if len(ccs) > 1:
        giant_comp = ccs[0]
        for cc in ccs[1:]:
            u,w = random.choice(giant_comp), random.choice(cc)
            G.add_edge(u,w)
            new_edges.add((u,w))
    return new_edges

def clean_c_data(G, c_data):
#this method serves no function other than trapping bugs: it removes data which should not be used in uncoarsening 
    aggregates    = c_data['aggregates']
    trapped_edges = c_data['trapped_edges']
    home_nodes    = c_data['home_nodes']
    merged_edges  = c_data['merged_edges']
    
    deleted_seeds   = [node for node in aggregates   if not G.has_node(node)]
    deleted_c_edges = [edge for edge in merged_edges if not G.has_edge(*edge)]

    for node in deleted_seeds:
        trapped_edges.pop(node)
        for guest in aggregates[node]:
            home_nodes.pop(guest)
        aggregates.pop(node)
    for edge in deleted_c_edges:
        merged_edges.pop(edge)

    return c_data

def compute_tpl_data(G, level, params):
#measures statistics of the topology of graph G
    tpl_data = {}
    tpl_data['enforce_connected'] = nx.is_connected(G)
    
    #estimates the probability of friending a node at distance d
    new_edge_horizon     = params.get('new_edge_horizon', estimate_horizon(G))  #no edges added to nodes beyond the horizon
    num_pairs_to_sample  = params.get('num_pairs_to_sample', 100)  #no edges added to nodes beyond the horizon
    locality_algorithm   = params.get('locality_algorithm', 'bypass_path_only')  #which method to use for computing locality
    triangle_distribution_limit = params.get('triangle_distribution_limit', None)  #truncation of the triangle distribution

    num_nodes_beyond_the_horizon = 0.
    overall_estimates    = np.zeros(new_edge_horizon+1)
    raw_mutual_nb_distribution = {}

    #for each node u in a sample, select one neighbor, and compute the distance up to H steps
    #    now see how many of the other neighbors of u have been reached.
    #    those not reached could be a kind of "beyond-the-horizon-edges", which are also possible
    for source in random.sample(G.nodes(), min(G.number_of_nodes(), num_pairs_to_sample)):
        source_degree = G.degree(source)
        if source_degree == 0:
            continue
        elif source_degree == 1:
            target_found = False
            #source_degree==1 is an important indicator that many edges are, in effect, chance edges
        else:
            target_nb = random.choice(G.neighbors(source))
            update_of_estimates, target_found, num_mutual_nbs = chart_paths(G, source, target_nb, new_edge_horizon, search_method=locality_algorithm, params=params)
            if num_mutual_nbs != None:
                if num_mutual_nbs in raw_mutual_nb_distribution:
                    raw_mutual_nb_distribution[num_mutual_nbs] += 1
                else:
                    raw_mutual_nb_distribution[num_mutual_nbs] = 1
        if target_found:
            overall_estimates += update_of_estimates
        else:
            num_nodes_beyond_the_horizon += 1
    try:
        locality_bias_correction = params['locality_bias_correction'][level]
    except:
        locality_bias_correction = 0.
    if locality_bias_correction > 0: #shift weight downward b/c this estimator under-rates correlations between neighbors
        overall_estimates[-1]          += locality_bias_correction     * num_nodes_beyond_the_horizon
        num_nodes_beyond_the_horizon   *= (1-locality_bias_correction)
        for dis in range(len(overall_estimates)-1, 2, -1):
            overall_estimates[dis-1] += locality_bias_correction     * overall_estimates[dis]
            overall_estimates[dis]   *= (1-locality_bias_correction)
    else:                            #shift weight upwards b/c this estimator over-rates correlations between neighbors
        for dis in range(len(overall_estimates)-1):
            overall_estimates[dis+1]   += -locality_bias_correction     * overall_estimates[dis]
            overall_estimates[dis]     *= (1+locality_bias_correction)
        num_nodes_beyond_the_horizon += -locality_bias_correction     * overall_estimates[-1]
        overall_estimates[-1]        *= (1+locality_bias_correction)

    accept_chance_edges = params.get('accept_chance_edges', 1.0)
    assert accept_chance_edges >= 0 and accept_chance_edges <= 1.0
    if sum(overall_estimates) > 0 or (num_nodes_beyond_the_horizon > 0 and accept_chance_edges > 0):
        if accept_chance_edges > 0:
            norm = accept_chance_edges*num_nodes_beyond_the_horizon + sum(overall_estimates)
            chance_edge_prob  = float(accept_chance_edges*num_nodes_beyond_the_horizon)/norm
        else:
            norm = sum(overall_estimates)
            chance_edge_prob  = 0.
        locality_acceptor = overall_estimates/norm
    else: #fallback
        locality_acceptor = [0., 0.] + [0.5/(2**d) for d in range(1, min(new_edge_horizon,G.number_of_nodes()-2))] 
        chance_edge_prob  = 0.
        if G.number_of_edges() > 10 and nx.density(G) > 0.2:
            print_warning(params, 'Warning: unable to estimate edge locality.')
            print_warning(params, 'Consider setting allow_chance_edges to positive values')
    assert locality_acceptor[0] == 0.
    assert locality_acceptor[1] == 0.
    if locality_bias_correction > 0 and locality_acceptor[2] > 0.8:
        print_warning(params, 'Warning: extreme locality at distance 2.  Might make it difficult to insert edges')

    #wishlist: need kernel density estimation for intermediate and low triangle densities
    if raw_mutual_nb_distribution != {}:
        if triangle_distribution_limit == None:
            triangle_distribution_limit = max(raw_mutual_nb_distribution.keys())
        mutual_nb_distribution = np.zeros(triangle_distribution_limit+2)  #the last element indicates the truncation
        for k in raw_mutual_nb_distribution:
            if k < triangle_distribution_limit:
                mutual_nb_distribution[k] = raw_mutual_nb_distribution[k]
            else:
                mutual_nb_distribution[-1] += 1.0
        mutual_nb_distribution /= mutual_nb_distribution.sum() + 1E-5
    else: #fallback
        mutual_nb_distribution = np.array([1.])
    
    tpl_data['locality_acceptor'] = locality_acceptor
    tpl_data['chance_edge_prob']  = chance_edge_prob 
    tpl_data['mutual_nb_distribution'] = mutual_nb_distribution

    return tpl_data

def do_coarsen(G, params):
    G_coarse = nx.empty_graph()
    aggregates = {} #nodes within new nodes.  seed->fine_nodes
    trapped_edges = {} #edges within new nodes.  seed->fine_edges
    home_nodes    = {} #node->seed
    merged_edges  = {} #edge->internal edges
    
    algorithm_for_coarsening      = params.get('algorithm_for_coarsening', seed_finder_weight_alg)
    seeds, home_nodes, aggregates = algorithm_for_coarsening(G, params)

    free_edges = set()        #edges not within any coarse node.  they will be retained in the coarse graph (many-to-one mapping)
    for seed in seeds:
        G_coarse.add_node(seed)
        G_coarse.node[seed]['weight'] = sum(G.node[nb].get('weight', 1.) for nb in aggregates[seed])

        trapped_edges[seed] = G.subgraph(aggregates[seed]).edges(data=False)

        for nb in aggregates[seed]:
            for nbnb in G.neighbors(nb):
                if nbnb in aggregates[seed] or (nbnb,nb) in free_edges:
                    continue
                free_edges.add((nb,nbnb))  

    for u,v in free_edges:
        s1 = home_nodes[u]
        s2 = home_nodes[v]
        uv_edge_wt = G.edge[u][v].get('weight', 1.0)
        if (s1,s2) in merged_edges:
            merged_edges[(s1,s2)].append((u,v))
            G_coarse.edge[s1][s2]['weight'] += uv_edge_wt
        elif (s2,s1) in merged_edges:
            merged_edges[(s2,s1)].append((u,v))
            G_coarse.edge[s2][s1]['weight'] += uv_edge_wt
        else:
            G_coarse.add_edge(s1,s2, weight=uv_edge_wt)
            merged_edges[(s1,s2)] = [(u,v)]
        assert (v,u) not in merged_edges[(s1,s2)]
    
    c_data = {'aggregates':aggregates, 'trapped_edges':trapped_edges, 'home_nodes':home_nodes, 'merged_edges':merged_edges}
    if 'do_coarsen_tester' in params:
        params['do_coarsen_tester'](G, G_coarse, c_data)

    return G_coarse, c_data


def do_uncoarsen(G_coarse, c_data, params):
    if isinstance(params.get('algorithm_for_uncoarsening', False), collections.Callable):
        return params['algorithm_for_uncoarsening'](G_coarse, c_data, params)

    aggregates    = c_data['aggregates']
    trapped_edges = c_data['trapped_edges']
    home_nodes    = c_data['home_nodes']
    merged_edges  = c_data['merged_edges']

    G_fine = nx.empty_graph()
    G_fine.add_nodes_from(home_nodes)

    for seed in trapped_edges:
        for u,v in trapped_edges[seed]:
            if u in G_fine and v in G_fine:
                G_fine.add_edge(u,v)
            #u or v must have been deleted
            
    for s1,s2 in G_coarse.edges_iter():
        if (s1,s2) in merged_edges:
            s1s2 = merged_edges[(s1,s2)]
        else:
            s1s2 = merged_edges[(s2,s1)]
            
        for u,v in s1s2:
            assert u in G_fine
            assert v in G_fine
            G_fine.add_edge(u,v)

    if 'do_uncoarsen_tester' in params:
        params['do_uncoarsen_tester'](G_coarse, G_fine, c_data)

    return G_fine




def edit_edges_sequential(G, edge_edit_rate, edge_growth_rate, tpl_data, params):
    #edit edges: first delete, then insert
    verbose = params.get('verbose', True)
    try:
        edit_rate = edge_edit_rate != [] and float(edge_edit_rate[0]) or 0.
        if edit_rate < 0. or edit_rate > 1.: raise
    except:
        print_warning(params, 'Bad or truncated edge edit rate information!  Defaulting to 0')
        edit_rate = 0.
    try:
        growth_rate = edge_growth_rate != [] and float(edge_growth_rate[0]) or 0.
    except:
        print_warning(params, 'Bad or truncated edge growth rate information!  Defaulting to 0')
        growth_rate = 0.
    if verbose:
        print('  Edge rates: edit %f, growth %f'%(edit_rate,growth_rate))
    if G.number_of_nodes() == 0:
        if verbose:
            print('Num nodes = 0 ... editing canceled')
        return G

    new_edge_horizon   = params.get('new_edge_horizon', estimate_horizon(G))  #no edges added to nodes beyond the horizon
    if new_edge_horizon in params and nx.density(G) > 0.2 and G.number_of_nodes() > 500:
        print_warning(params, 'Warning: using a large horizon (%d) on a large graph might use a lot of time'%new_edge_horizon)

    if 'enforce_connected' in params:
        enforce_connected = params['enforce_connected'] 
    else:
        enforce_connected = tpl_data['enforce_connected']
    dont_cutoff_leafs = params.get('dont_cutoff_leafs', False)
    #do we allow leafs to be cut off completely?
    #   this option should be used sparingly, as it disrupts deferential retachment and decreases clustering

    all_nodes = G.nodes() 
    added_edges_set = set()
    deled_edges_set = set()
    target_edges_to_delete = npr.binomial(max(G.number_of_edges(), 1), edit_rate)
    target_edges_to_add    = npr.binomial(max(G.number_of_edges(), 1), edit_rate)  #should be here, since NumEdges will change
    if growth_rate > 0:
        target_edges_to_add    += int(round(G.number_of_edges() * growth_rate))
    else:
        target_edges_to_delete += int(round(G.number_of_edges() * (-growth_rate)))

    deprived_nodes = [] #list of nodes that lost edges, including repetitions
    deferential_detachment_factor = params.get('deferential_detachment_factor', 0.0)  
    avg_degree = np.average(list(nx.degree(G).values()))  #inexact deferential detachment, but with a much higher sampling efficiency
    num_deletion_trials           = params.get('num_deletion_trials', int(round(avg_degree**2)) )

    G_adj = G.adj
    G_degree    = lambda u: G_adj[u].__len__()
    G_neighbors = lambda u: list(G_adj[u].keys())
    for trial_num in range(max(20, num_deletion_trials*target_edges_to_delete)):
        if len(deled_edges_set) == target_edges_to_delete:
            break
        u = random.choice(all_nodes)
        degree_of_u = G_degree(u)
        if degree_of_u == 0: #will take care of this later
            continue 
        w = random.choice(G_neighbors(u))
        degree_of_w = G_degree(w)
        #perhaps a slight improvement is to multiply not by avg_degree but by avg_nb_degree
        if npr.rand()*deferential_detachment_factor > avg_degree/float(degree_of_u*degree_of_w):
            continue
        if dont_cutoff_leafs and (degree_of_u == 1 or degree_of_w == 1):
            continue
        if strong_clustering_structure(G, u, w, params):
            continue
        G.remove_edge(u,w)
        deled_edges_set.add((u,w))
        deprived_nodes += [u,w]

    if enforce_connected: 
        new_edges = check_and_fix_connectivity(G, params)
        added_edges_set.update(new_edges)
    
    num_remaining_edges_to_add = target_edges_to_add - len(added_edges_set)
    edge_welfare_fraction = params.get('edge_welfare_fraction', 0.0)  
    long_bridging         = params.get('long_bridging', False)  
    #whether it should try to build edges to nodes which lost them;  not supported for all edges or for edges lost during node deletion

    for trial_num in range(max(20, 3*target_edges_to_add)): 
        if num_remaining_edges_to_add <= 0:  #we might overshoot, hence <= 0 not ==
            break
        if npr.rand() > edge_welfare_fraction or len(deprived_nodes) == 0:
            head = random.choice(all_nodes)
        else:
            head = random.choice(deprived_nodes)
        if G_degree(head) == 0 and tpl_data['chance_edge_prob'] == 0.0:
            continue 
        distance_from_head = graphutils.bfs_distance_with_horizon(G, source=head, horizon=new_edge_horizon, blocked_node=None)
        tail = find_node_to_friend(G=G, head=head, distance_from_head=distance_from_head, tpl_data=tpl_data, params=params, existing_nbs=G_neighbors(head))
        if tail == None or tail == head or G.has_edge(head,tail):
            continue
        head_to_tail_distance = distance_from_head.get(tail, -1)  #-1 is the case of infinity
        if long_bridging and head_to_tail_distance>1:
            back_node = head
            for dummy in range(head_to_tail_distance-1):
                new_node = new_node_label(G) 
                G.add_edge(back_node,new_node)  #this also adds a node
                added_edges_set.add((back_node,new_node))
                back_node = new_node
            G.add_edge(back_node,tail)
            added_edges_set.add((new_node,tail))
            num_remaining_edges_to_add -= head_to_tail_distance
        else:
            G.add_edge(head,tail)
            added_edges_set.add((head,tail))
            num_remaining_edges_to_add -= 1
        
    num_edges_added   = len(added_edges_set)
    num_edges_deleted = len(deled_edges_set)
    #print num_edges_added, num_edges_deleted, G.number_of_edges()
    if num_edges_added > 20 and (num_edges_added-target_edges_to_add)/float(num_edges_added)> 0.2:
        print_warning(params, 'Warning: Excessive number of edges were added. Is the graph treelike (low AvgDegree and connected)? AvgDegree=%.1f.'%np.average(list(nx.degree(G).values())))
        #this might be caused by node edits.  in that case, try minorizing_node_deletion
    if num_edges_added > 20 and (target_edges_to_add-num_edges_added)/float(num_edges_added)> 0.2:
        print_warning(params, 'Warning: Excessive number of edges failed to add. Is the graph too dense? Density=%.2f'%nx.density(G))
    if num_edges_deleted > 20 and abs(target_edges_to_delete-num_edges_deleted)/float(num_edges_deleted)> 0.2:
        print_warning(params, 'Warning: Excessive number of edges were deleted. Is the graph too dense? Density=%.2f'%nx.density(G))
    if verbose:
        print('\tadded edges: %d, deleted edges: %d'%(num_edges_added,num_edges_deleted))

    if 'edit_edges_tester' in params:
        params['edit_edges_tester'](G, added_edges_set, deled_edges_set, tpl_data)


    return G



def edit_nodes_sequential(G, node_edit_rate, node_growth_rate, tpl_data, params):
    verbose = params.get('verbose', True)
    if verbose:
        print('nn: %d'%G.number_of_nodes())
    try:
        edit_rate = node_edit_rate != [] and float(node_edit_rate[0]) or 0.
        if edit_rate < 0. or edit_rate > 1.: raise
    except:
        print_warning(params, 'Bad or truncated node edit rate information!  Defaulting to 0')
        edit_rate = 0.
    try:
        growth_rate = node_growth_rate != [] and float(node_growth_rate[0]) or 0.
    except:
        print_warning(params, 'Bad or truncated node growth rate information!  Defaulting to 0')
        growth_rate = 0.
    if verbose:
        print('  Node rates: edit %f, growth %f'%(edit_rate,growth_rate))
    if G.number_of_nodes() == 0:
        if verbose:
            print('Num nodes = 0 ... editing canceled')
        return G

    new_edge_horizon   = params.get('new_edge_horizon', estimate_horizon(G))  #no edges added to nodes beyond the horizon
    if new_edge_horizon in params and new_edge_horizon > 5 and nx.density(G) > 0.2 and G.number_of_nodes() > 500:
        print_warning(params, 'Warning: using a large horizon (%d) on a large graph might use a lot of time'%new_edge_horizon)
    
    num_deleted_nodes = npr.binomial(G.number_of_nodes(), edit_rate)
    num_added_nodes   = npr.binomial(G.number_of_nodes(), edit_rate) 
    if growth_rate > 0:
        num_added_nodes   += int(round(G.number_of_nodes() * growth_rate))
    else:
        num_deleted_nodes += int(round(G.number_of_nodes() * (-growth_rate)))
    if num_deleted_nodes > G.number_of_nodes():
        print_warning(params, 'Warning: excess negative growth rate. Deletion of nodes will destroy all the nodes of the graph.  Editing aborted at this level.')
        return G

    G_adj = G.adj
    G_degree    = lambda u: G_adj[u].__len__()
    G_neighbors = lambda u: list(G_adj[u].keys())
    #we cache edges-to-add to avoid densifying the graph during the loop
    original_nodes = G.nodes()
    added_node_info = {} 
    for i in range(num_added_nodes):
        source_node = random.choice(original_nodes)
        new_node    = new_node_label(G)
        added_node_info[new_node] = G_degree(source_node)
        #G.node[new_node]['resampling_source'] = source_node

    num_edges_added = 0
    num_edges_deleted = 0
    failed_searches = 0
    for new_node in added_node_info:
        G.add_node(new_node)
        num_remaining_nbs_to_add = added_node_info[new_node]
        if num_remaining_nbs_to_add == 0:
            continue
        #uncomment below and use 'enforce_connected':False to see the aggregates.  WARNING: comment back when done
        #continue

        anchor_node = random.choice(original_nodes)
        G.add_edge(new_node, anchor_node)
        num_edges_added += 1
        num_remaining_nbs_to_add -= 1
        for trial_num in range(max(40, 3*num_remaining_nbs_to_add)):
            if num_remaining_nbs_to_add == 0:
                break
            distance_new_node = graphutils.bfs_distance_with_horizon(G, source=new_node, horizon=new_edge_horizon)
            v = find_node_to_friend(G=G, head=new_node, distance_from_head=distance_new_node, tpl_data=tpl_data, params=params, existing_nbs=G_neighbors(new_node), offset=0)
            if v == None or v == new_node or v in G_neighbors(new_node):
                continue
            G.add_edge(new_node, v)
            num_edges_added += 1
            num_remaining_nbs_to_add -= 1
        if num_remaining_nbs_to_add > 0:
            failed_searches += 1
    added_nodes_set = set(added_node_info.keys())
        
    deled_nodes_set = set()
    minorizing_node_deletion = params.get('minorizing_node_deletion', False)
    for u in random.sample(G.nodes(), num_deleted_nodes):
        num_edges_deleted += G_degree(u)
        if minorizing_node_deletion: #connect the neighbors into a tree
            nbrs = G_neighbors(u)
            random.shuffle(nbrs)
            for nb_idx,nb in enumerate(nbrs[:-1]):  
                new_edge = (nb, nbrs[nb_idx+1])
                if not G.has_edge(*new_edge):
                    #assert new_edge[0] != new_edge[1]
                    G.add_edge(*new_edge)
                    num_edges_added += 1
        G.remove_node(u)
        deled_nodes_set.add(u)

    if len(original_nodes) !=0 and len(original_nodes)>10 and (float(failed_searches)/len(original_nodes)) > .2:
        print_warning(params, 'Warning: > 20%% of searches failed when attempting to insert edges. Is the graph too dense? Density=%.2f'%nx.density(G))
    if verbose:
        print('\tadded nodes: %d, deleted nodes: %d'%(num_added_nodes,num_deleted_nodes))
        print('\tadded edges: %d, deleted edges: %d'%(num_edges_added,num_edges_deleted))

    if 'edit_nodes_tester' in params:
        params['edit_nodes_tester'](G, added_nodes_set, deled_nodes_set, tpl_data)

    return G

def estimate_horizon(G):
    #how much should we branch if we don't want to visit more than 100 nodes during BFS?
    density = nx.density(G)

    if density == 0 or G.number_of_nodes() < 3:
        return 4

    new_nbs_per_level = (G.number_of_nodes() - 1) * density  
    if int(round(new_nbs_per_level)) == 1:
        return 4
    #assuming a random graph, we reach 100 nodes in the final level if the horizon satisfies: 
    #100 = new_nbs_per_level ** horizon 
    my_estimate = int(round(np.log(100)/np.log(new_nbs_per_level)))

    return max(2, my_estimate)

def find_node_to_friend(G, head, distance_from_head, tpl_data, params, existing_nbs=None, offset=0):
#find nodes to friend
    locality_acceptor      = tpl_data['locality_acceptor'] 
    mutual_nb_distribution = tpl_data['mutual_nb_distribution']
    #implicit: chance_edge_prob  = tpl_data['chance_edge_prob'] 

    all_nodes = None
    tail = None
    if existing_nbs == None:
        existing_nbs = set()
    else:
        existing_nbs = set(existing_nbs)
    lookup_by_dis = {}
    G_adj = G.adj
    G_neighbors = lambda u: list(G_adj[u].keys())
    for candidate in distance_from_head:
        dis = distance_from_head[candidate]+offset
        if dis not in lookup_by_dis:
            lookup_by_dis[dis] = []
        if candidate in existing_nbs:
            continue
        lookup_by_dis[dis].append(candidate)

    num_insertion_trials = params.get('num_insertion_trials', 10 )
    for trial in range(num_insertion_trials):
        #sample from np.random.multinomial()
        toss = npr.rand() 
        for dis, prob in enumerate(locality_acceptor):
            toss -= prob
            if toss < 0:
                break
        if toss < 0:
            if dis not in lookup_by_dis:
                continue  #might occur legitimately when we have offset=1
            candidates = lookup_by_dis[dis]
            #wishlist: we could build lookup_by_dis lazily, only if we know dis
            if len(candidates) == 0:
                continue
            if dis != 2 or not params.get('fine_clustering', False):
                candidate = random_one(candidates)
                tail = candidate
            else:
                scores   = []
                head_nbs = set(G_neighbors(head))
                for candidate in candidates:
                    num_mutual_nbs = len(head_nbs.intersection(G_neighbors(candidate)))
                    if num_mutual_nbs < len(mutual_nb_distribution)-1: 
                        cand_score = mutual_nb_distribution[num_mutual_nbs]
                    else:
                        cand_score = mutual_nb_distribution[-1]
                    scores.append(cand_score)
                scores_sum = sum(scores)
                if scores_sum > 0:
                    scores = np.array(scores)/scores_sum
                    candidate_idx = npr.multinomial(1, scores).nonzero()[0][0]
                    tail = candidates[candidate_idx]
                else:
                    tail = random.choice(candidates)
                #print scores[candidate_idx]
                #print len(head_nbs.intersection(G_neighbors(candidate)))
                #print mutual_nb_distribution
                #print
        else:
            for i in range(G.number_of_nodes()):
                if all_nodes == None:
                    all_nodes = G.nodes()
                candidate = random.choice(all_nodes)
                if (candidate not in distance_from_head) and candidate != head and (candidate not in existing_nbs):
                    tail = candidate
                    break
        if tail != None:
            break
    return tail


def interpolate_edges(G, c_data, model_map, fine_model_map, params):    
    aggregates    = c_data['aggregates']
    merged_edges  = c_data['merged_edges']
    deep_copying  = params.get('deep_copying', True)
    if not deep_copying: assert len(model_map) == 0
    
    authentic_edges = list(merged_edges.items())
    edited_edges = []
    for (s1,s2) in G.edges_iter(): #wishlist: faster loop?
        if ((s1,s2) not in merged_edges) and ((s2,s1) not in merged_edges):
            edited_edges.append((s1,s2))
    for (s1,s2) in edited_edges:
        trapped_in_s1 = aggregates[s1]
        trapped_in_s2 = aggregates[s2]
        
        model_aggregate_1 = model_map.get(s1, None)
        model_aggregate_2 = model_map.get(s2, None)
        new_pairs = set()
        if deep_copying and G.has_edge(model_aggregate_1, model_aggregate_2):
            if (model_aggregate_1, model_aggregate_2) in merged_edges:
                merged_model_edge_contents = merged_edges[(model_aggregate_1, model_aggregate_2)]
            else:
                merged_model_edge_contents = merged_edges[(model_aggregate_2, model_aggregate_1)]
            reversed_map = {}
            for trapped_node in trapped_in_s1 + trapped_in_s2:
                reversed_map[fine_model_map[trapped_node]] = trapped_node
            for mA, mB in merged_model_edge_contents:
                u = reversed_map[mA]
                v = reversed_map[mB]
                new_pairs.add((u,v))
        else:
            if authentic_edges != []:
                random_model_edge, random_model_edge_contents = random.choice(authentic_edges)
                num_target_edges = len(random_model_edge_contents)
            else:
                num_target_edges = 1
            num_failures = 0
            while len(new_pairs) < num_target_edges and num_failures < max(10, 3*num_target_edges):
                u = random.choice(trapped_in_s1)
                v = random.choice(trapped_in_s2)
                if ((u,v) not in new_pairs) and ((v,u) not in new_pairs):
                    new_pairs.add((u,v))
                else:
                    num_failures += 1
        merged_edges[(s1,s2)] = [(u,v) for u,v in new_pairs]
        #might optionally pass an attribute 'new' in the edges
    return c_data


def interpolate_nodes(G, c_data, model_map, params):    
#data in aggregates from from the graph before editing.  we need to augment it with interpolation of new nodes.
#interpolation is of two kinds: (base case) G just now received a new node -> we resample.  (deep copying) a new node appeared at earlier levels, leading to interpolation now.
#the distinction between model_map and fine_model_map is this: the former are clones, while the latter are children of clone nodes
    aggregates    = c_data['aggregates']
    trapped_edges = c_data['trapped_edges']
    home_nodes    = c_data['home_nodes']
    merged_edges  = c_data['merged_edges']
    if not params.get('deep_copying', True):
        assert len(model_map) == 0

    #this_level_node A -> node in G_i original of which is the model of A 
    fine_model_map = {}
   
    authentic_nodes = list(aggregates.keys())
    assert authentic_nodes != []
    edited_nodes = [node for node in G if node not in aggregates]
    num_new_nodes = 0
    num_new_edges = 0
    for node in edited_nodes:
        source_aggregate = model_map.get(node, random.choice(authentic_nodes))
        sources_edges    = trapped_edges[source_aggregate]
        sources_nodes    = aggregates[source_aggregate]

        renamed_nodes = {}
        for node_hosted_by_source in sources_nodes:
            new_hosted_node                        = new_node_label(home_nodes)
            renamed_nodes[node_hosted_by_source]   = new_hosted_node
            fine_model_map[new_hosted_node]        = node_hosted_by_source
            num_new_nodes += 1

        my_trapped_nodes = list(renamed_nodes.values())

        my_trapped_edges = []
        for edge_hosted_by_source in sources_edges:
            my_trapped_edges.append( (renamed_nodes[edge_hosted_by_source[0]], renamed_nodes[edge_hosted_by_source[1]]) )
            num_new_edges += 1

        aggregates[node]    = my_trapped_nodes
        trapped_edges[node] = my_trapped_edges
        for new_hosted_node in my_trapped_nodes:
            home_nodes[new_hosted_node] = node
        #might optionally pass an attribute 'new' in the node
        #num_added_nodes += len(my_trapped_nodes)
    if params.get('verbose', True):
        print('  from new aggregates: %d nodes, %d edges'%(num_new_nodes,num_new_edges))
    #print 'added: %d'%num_added_nodes

    if params.get('deep_copying', True):
        return c_data, fine_model_map
    else:
        return c_data, {}

def flush_graph(G):            
#the algorithm relabels the nodes of the graph at random 
#-> data on aggregates becomes obsolete and is so never used
    node_map = {} 
    for node in G:
        while True:
            new_name = new_node_label(G)
            if (new_name not in G) and (new_name not in node_map):
                break
        node_map[node] = new_name
    G = nx.relabel_nodes(G, node_map, copy=True)
    return G

def generate_graph(original, params=None):
    if params == None:
        params = {}
        print_warning(params, 'WARNING: empty parameter input. Running with default parameters.')

    if params.get('algorithm', False):
        params2 = params.copy()
        alg_info = params2.pop('algorithm')
        if isinstance(alg_info, collections.Callable):
            alg_method = alg_info
        elif type(alg_info) is str:
            alg_method = eval(alg_method)
        elif (type(alg_info) is list) or (type(alg_info) is tuple):
            alg_method           = alg_info[0]
            params2['algorithm'] = alg_info[1]
        else:
            raise ValueError('algorithm parameter should be either callable, the name of a function, or (func,(nested_algorithm))')
        return alg_method(original=original, params=params2)

    simpletesters.validate_params(params)

    node_edit_rate   = params.get('node_edit_rate', [])
    edge_edit_rate   = params.get('edge_edit_rate', [])
    node_growth_rate = params.get('node_growth_rate', [])
    edge_growth_rate = params.get('edge_growth_rate', [])
 
    #we might want to convert all nodes to integers for performance reasons

    G = original.copy()
    if params.get('verbose', True):
        print('Checking original graph ...')
        simpletesters.sanity_test(G, params)
    start_time = time.time()

    replica, model_map = revise_graph(G=G, level=0,
                                node_edit_rate=node_edit_rate, 
                                node_growth_rate=node_growth_rate, 
                                edge_edit_rate=edge_edit_rate, 
                                edge_growth_rate=edge_growth_rate, 
                                params=params)
    if params.get('verbose', True):
        print('replica is finished. nn: %d.  time: %.2f sec.'%(replica.number_of_nodes(), time.time()-start_time))
        print()
    
    replica      = resample_attributes(G, replica, model_map, params)
    replica.name = getattr(original, 'name', 'graph') + '_replica_' + timeNow()

    simpletesters.sanity_test(replica, params)
    del G
    gc.collect()

    return replica


def musketeer_on_subgraphs(original, params=None):
    components = nx.connected_component_subgraphs(original)
    merged_G   = nx.Graph()

    component_is_edited = params.get('component_is_edited', [True]*len(components))

    for G_num, G in enumerate(components):
        if component_is_edited[G_num]:
            replica = generate_graph(original=G, params=params)
        else:
            replica = G
        
        merged_G = nx.union(merged_G, replica)

    merged_G.name = getattr(original, 'name', 'graph') + '_replica_' + timeNow()
    
    return merged_G

def musketeer_snapshots(original, params=None):
#applies replication sequentially, generating snapshots of the original.  
#returns the final snapshot.  
#snapshots (0 to last) are held in the .snapshot attribute
    graphs = [original]

    num_snapshots = params['num_snapshots']

    for graph_num in range(num_snapshots):
        G = graphs[-1]
        replica = generate_graph(original=G, params=params)
        replica.name = 'snapshot_%d'%graph_num
        graphs.append(replica)

    if params.get('verbose', True):
        print('Snapshots complete.')
        print()

    replica.snapshots = graphs
    
    return replica


def musketeer_iterated_cycle(original, params=None):
#applies replication sequentially, and returns the final snapshot.  
    num_cycles = params['num_v_cycles']

    params2 = params.copy()
    params2['edge_edit_rate']   = [r/float(num_cycles) for r in params2.get('edge_edit_rate', [])]
    params2['edge_growth_rate'] = [r/float(num_cycles) for r in params2.get('edge_growth_rate', [])]
    params2['node_edit_rate']   = [r/float(num_cycles) for r in params2.get('node_edit_rate', [])]
    params2['node_growth_rate'] = [r/float(num_cycles) for r in params2.get('node_growth_rate', [])]
    
    replica = original
    for graph_num in range(num_cycles):
        replica = generate_graph(original=replica, params=params2)
        replica.name = getattr(original, 'name', 'graph') + '_replica_w%d_'%graph_num + timeNow()

    return replica 


def new_node_label(G):
#G is either a graph or a dict/list of existing labels
    num_trials = 100
    label = None
    for t in range(num_trials):
        label = npr.randint(max_int)
        if label not in G:
            break
    if label == None:
        raise Exception('Could not find a unique label for a newly-added node')
    return label

def print_warning(params, str):
    if not params.get('suppress_warnings', False): 
        print(str)

def resample_attributes(G, replica, model_map, params):
#inserts attributes to new nodes and edges by COPYING data from existing nodes and edges
#note that to save the data, the replica should be saved as particular formats (e.g. gml, dot, edgelist[for edges only])
    maintain_node_attributes = params.get('maintain_node_attributes', False)
    maintain_edge_attributes = params.get('maintain_edge_attributes', False)
    deep_copying             = params.get('deep_copying', True)

    original_nodes = G.nodes()
    G_adj = G.adj
    G_degree    = lambda u: G_adj[u].__len__()
    G_neighbors = lambda u: list(G_adj[u].keys())
    if maintain_node_attributes:
        if params.get('verbose', True):
            print('Updating node attributes ...')
        for node in replica:
            if replica.node[node] != {}:
                continue
            elif node in G:
                replica.node[node] = G.node[node].copy()
            elif deep_copying:
                model_node = model_map.get(node, random.choice(original_nodes))
                replica.node[node] = G.node[model_node].copy()
            else:
                model_node = random.choice(original_nodes)
                replica.node[node] = G.node[model_node].copy()

    if maintain_edge_attributes and G.number_of_edges() > 0:
        if params.get('verbose', True):
            print('Updating edge attributes ...')
        for edge in replica.edges_iter():
            if replica.get_edge_data(*edge) != {}:
                continue
            elif G.has_edge(*edge):
                edata = G.get_edge_data(*edge).copy()
                replica.edge[edge[0]][edge[1]] = edata
                replica.edge[edge[1]][edge[0]] = edata
            else:
                modelA = model_map.get(edge[0], None)
                modelB = model_map.get(edge[1], None)
                if deep_copying and G.has_edge(modelA, modelB):
                    nodeA, nodeB = modelA, modelB
                else:
                    for trial in range(G.number_of_edges()):
                        nodeA = random.choice(original_nodes)
                        if G_degree(nodeA) == 0:
                            continue
                        nodeB = random.choice(G_neighbors(nodeA))
                        break

                edata = G.get_edge_data(nodeA, nodeB).copy()
                replica.edge[edge[0]][edge[1]] = edata
                replica.edge[edge[1]][edge[0]] = edata

    return replica

def revise_graph(G, level, node_edit_rate, node_growth_rate, edge_edit_rate, edge_growth_rate, params):
#revises a graph at a particular resolution and deeper
    no_more_work   = node_edit_rate == [] and edge_edit_rate == [] and node_growth_rate == [] and edge_growth_rate == []
    excess_density = G.number_of_edges() == 0 or nx.density(G) > params.get('coarsening_density_limit', 0.9)
    if no_more_work or excess_density:
        if no_more_work and params.get('verbose', True): 
            print('No changes at this level or coarser. nn: %d'%G.number_of_nodes())
        if excess_density:
            if params.get('verbose', True): 
                print('Coarsening stopped due to excess density')
                if not no_more_work:
                    print('Editing at deeper levels is impossible.  CHANGE editing parameters.')
        if params.get('verbose', True): 
            print('-----------------')
        if params.get('memoriless_interpolation', False): 
            G = flush_graph(G)
        G.coarser_graph = None
        return G, {}
    
    #recursive operation on the coarse network.  model_map: nodeB ->nodeA of G which is model of nodeB
    coarse, c_data         = do_coarsen(G, params)
    coarse, model_map      = revise_graph(coarse, level+1, node_edit_rate[1:], node_growth_rate[1:], edge_edit_rate[1:], edge_growth_rate[1:], params)
    c_data, fine_model_map = interpolate_nodes(coarse, c_data, model_map, params=params)  #must precede edges, since otherwise cannot find an attachment for some edges
    c_data                 = interpolate_edges(coarse, c_data, model_map, fine_model_map, params=params)
    clean_c_data(coarse, c_data)
    G_prime = do_uncoarsen(coarse, c_data, params)
    
    #editing of the current level
    tpl_data = compute_tpl_data(G, level, params)

    G_prime  = edit_nodes_sequential(G_prime, node_edit_rate, node_growth_rate, tpl_data, params)
    G_prime  = edit_edges_sequential(G_prime, edge_edit_rate, edge_growth_rate, tpl_data, params)
    
    if 'revise_graph_tester' in params:
        params['revise_graph_tester'](G, G_prime, c_data)
    if params.get('retain_intermediates', False): 
        G_prime.model_graph, G_prime.coarser_graph = G, coarse
    else:
        G_prime.model_graph, G_prime.coarser_graph = None, None

    return G_prime, fine_model_map


def seed_finder_weight_alg(G, params):
    '''
    coarsening inspired by Safro et al. 
        Graph minimum linear arrangement by multilevel
        weighted edge contractions
    1. nodes have weights that represent nodes previously included in them
    2a. nodes with high expected weight (>= 2 of average) automatically become seeds
        expected weight = self weight + total weight of neighbors
    2b. we go over the remaining nodes in order of decreasing future weight
        and make them seeds when (sum of neighbor expect_wts who are seeds) / (sum over all nbs) < 0.4  
    3. the remaining nodes are assigned to a neighboring seeds that has the greatest wt
    '''

    seeds         = set()
    aggregates    = {}    #nodes within new nodes.  seed->fine_nodes
    home_nodes    = {}    #the reverse map: node->seed
    free_nodes    = set(G.nodes())
    seed_threshold_1 = params.get('seed_threshold_1', 2.0)
    seed_threshold_2 = params.get('seed_threshold_2', 0.4)

    G_adj = G.adj
    G_degree    = lambda u: G_adj[u].__len__()
    G_neighbors = lambda u: list(G_adj[u].keys())
    expected_wt_dict = {}
    for node in G:
        if G_degree(node) == 0: #singletons mess up seed finding by pulling down the average wt, but then become seeds themselves
            seeds.add(node)
        else:
            expected_wt_dict[node] = G.node[node].get('weight', 1.) + sum(G.node[nb].get('weight', 1.) for nb in G_neighbors(node))

    if G.number_of_edges() > 0:
        mean_exp_wt = np.average(list(expected_wt_dict.values())) 
    else:
        mean_exp_wt = 0.

    expected_wt = list(expected_wt_dict.items())
    expected_wt.sort(key=lambda x:x[1], reverse=True)
    for node,wt in expected_wt:
        if wt > seed_threshold_1*mean_exp_wt:
            seeds.add(node)
        elif sum(expected_wt_dict[nb] for nb in seeds.intersection(G_neighbors(node)))/sum(expected_wt_dict[nb] for nb in G_neighbors(node)) < seed_threshold_2:
            seeds.add(node)
    free_nodes.symmetric_difference_update(seeds)
    for seed in seeds:
        aggregates[seed] = [seed]
        home_nodes[seed] = seed

    for node in free_nodes:
        nearby_seeds = seeds.intersection(G_neighbors(node))
        seed_attraction = [(seed,G.edge[node][seed].get('weight', 1.0)) for seed in nearby_seeds] 
        my_aggregate = max(seed_attraction, key=lambda x:x[1])[0]
        aggregates[my_aggregate].append(node)
        home_nodes[node] = my_aggregate

    if params.get('verbose', True):
        print('nn: %d (seeds: %d)'%(G.number_of_nodes(),len(seeds)))
    return list(seeds), home_nodes, aggregates

def strong_clustering_structure(G, u, w, params):
#detects strong clustering around (u,w) and decides whether to skip deletion of the edge
    if not params.get('preserve_clustering_on_deletion', True):
        return False

    num_mutual_nbs = len(set(G.neighbors(u)).intersection(G.neighbors(w)))
    if num_mutual_nbs == 0:
        return False
    elif npr.rand() * num_mutual_nbs < 0.5:
        return False
    else:
        return True


