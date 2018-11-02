# -*- coding: utf-8 -*-
'''
Multiscale Entropic Network Generator (MUSKETEER)

Copyright (c) 2011-2014 by Alexander Gutfraind and Ilya Safro. 
All rights reserved.

Use and redistribution of this file is governed by the license terms in
the LICENSE file found in the project's top-level directory.


Benchmarks of Quality and performance - internal 


TODO:
    1. alg-alg benchmarks:
        compare two algorithms on a set of graphs.  
        a. attempt to show one is better on average using t-tests
        b. identify areas where one algorithm is consistently better, using t-tests

    2. target graph benchmarks:
        report how well a particular set of graphs is replicated by a given algorithm

    3. feature benchmark:
        a. one-step benchmark
        report how well a particular feature (e.g. clustering) is replicated by a given algorithm

        b. two-step benchmark:
        report how well a particular algorithm is able to retain the original features over 2 replications

        c. T-step benchmark
        report on the divergence from the original graph when replication is iterated

DONE    d. bias benchmark
        does the algorithm have a propensity to change a given feature in a particular direction? eg. increase density

    4. benchmark suites
        a set of tests the generates figures in the paper for PNAS

    5. noise-spectrum benchmarks
        show how changing the color of the noise changes the graph

    6. ensemble_validation_benchmark
        plot the properties of the replicas against properties of originals and against just randomly-noised graphs

    7. speed benchamarks
        running time

'''

import os, subprocess
import time
import numpy as np
import numpy.random as npr
import random, sys
import networkx as nx
import matplotlib
matplotlib.use('PS')
#matplotlib.use('PDF')
import matplotlib.pylab as pylab
#import pylab
import pdb, traceback
import pickle
import scipy.stats

#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#matplotlib.rc('font',**{'family':'serif','serif':['Palatino']})
#rc('font',**{'family':'serif','serif':['Times New Roman']})
try:
    os.system('latex -version > /dev/null')
    matplotlib.rc('text', usetex=True)
except:
    matplotlib.rc('text', usetex=False)

import algorithms
import alternatives
import community
import graphutils
import simpletesters, testers

np.seterr(all='raise')

timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())

if not os.path.exists('output'):
    os.makedir('output')

def alg_alg(alg1, alg2, params={}, graphs={}, alg_params={}, alg2_params={}):
    '''
    compare two algorithms
    alg1 and alg2 should be methods that accept a graph and alg_params;  
    alg2_params are the same as alg_params by default

    TODO: allow alg2 to be just random edge noise

    graphs[graph_name] = {'G':G, 'options':Not_implemented}
    '''
    if graphs == {}:
        graphs = next(default_graphs())
    if alg_params == {}:
        alg_params['verbose'] = False
        alg_params['error_rate'] = [0.05/(1.+i) for i in range(100)]
    if alg2_params == {}:
        alg2_params = alg_params

    errors_by_graph = []
    errors_by_type  = []
    errors_by_metric= []
    for graph_name in graphs:
        G          = graphs[graph_name]['G']
        graph_type = graphs[graph_name]['graph_type']

        replica1 = alg1(original=G, params=alg_params)
        replica2 = alg2(original=G, params=alg2_params)

        error_data1 = compare_nets(old_G=G, new_G=replica1, params={'verbose':False})[0]
        error_data2 = compare_nets(old_G=G, new_G=replica2, params={'verbose':False})[0]
        errors_by_graph.append((graph_name,
                                np.average([abs(err) for err in list(error_data1.values())]),
                                np.average([abs(err) for err in list(error_data2.values())])))

    print('Mean errors by graph (absolute)')
    print('Graph\t\t\t\tAlg1_Err Alg2_Err')
    for graph_name, err1, err2 in errors_by_graph:
        print('%s\t%.2f\t%.2f'%(graph_name.center(25),err1,err2))
    print()

    print('TODO: report mean errors by graph type')
    print('GraphType\t\tAlg1_Error\t\tAlg2_Error')
    print()

    print('TODO: report errors by metric')
    print('Metric\t\tAlg1_Error\t\tAlg2_Error')
    print()

def alternatives_metric(G=None, replicas=None, metric=(community.louvain_modularity, 'modularity'), seed=None, figpath=None, title_infix=''):
    vals_of_replicas = {}
    num_of_algs = len(replicas)
    alg_names   = list(replicas.keys())
    for alg_name in alg_names:
        print(alg_name)
        vals_of_replicas[alg_name] = []
        for num_replica, replica in enumerate(replicas[alg_name]):
            vals_of_replicas[alg_name].append(metric['function'](replica))
            #if num_replica > 10:
            #    break
            sys.stdout.write('.')
            sys.stdout.flush()
        sys.stdout.write('; ')
    print()
    val_of_graph = metric['function'](G)
    print('Val of graph: ')
    print(val_of_graph)
    
    med_vals = [np.median(vals_of_replicas[alg_name]) for alg_name in alg_names]
    print('Medians:')
    print(med_vals)

    normed_replica_vals = []
    for alg_idx, alg_name in enumerate(alg_names):
        noravg_vals = np.array(vals_of_replicas[alg_name])/((1E-20) + val_of_graph)
        normed_replica_vals.append(noravg_vals)

    fig = pylab.figure()
    pylab.hold(True)
    pylab.plot(num_of_algs*[1.], list(range(num_of_algs)), 'o', color='k', linewidth=2., label=G.name)
    
    pylab.yticks(list(range(num_of_algs)), [alg_name.replace(' ', '\n') for alg_name in alg_names], rotation=0)
    pylab.xlabel(metric['name'] + ' (normalized by empirical value)', rotation=0, fontsize='15')
    #pylab.xlabel('Normalized value', rotation=0, fontsize='20')#, x=0.1)
    pylab.boxplot(np.array(normed_replica_vals).transpose(), positions=list(range(num_of_algs)), vert=0, patch_artist=True)
    max_axis = 2
    pylab.xlim(-0.02,max_axis)
    fig.subplots_adjust(left=0.17, right=0.95)
    '''
    pylab.text(x=1.75, y=num_of_algs-0.4, s='Median of\nreplicas')
    for alg_idx, alg_name in enumerate(alg_names):
        normed_val = med_vals[alg_idx]/((1E-20) + val_of_graph)
        if abs(normed_val) < 100:
            val_str=r'$%.2f$'%normed_val
        elif normed_val < -1000:
            val_str='un-\ndefined'
        else:
            val_str=r'$\gg0$'
        pylab.text(x=max_axis+0.02, y=alg_idx, s=val_str)
    '''
    pylab.hold(False)

    if figpath == None:
        figpath = 'output/alternatives_metric_'+metric['name']+'_'+G.name+'__'+timeNow()
        figpath = clean_path(figpath)
    save_figure_helper(figpath)


def alternatives_outbreaks(G=None, seed=None, figpath=None, params=None, dynamics_params=None, alternative_algs=None, num_replicas = 150, outbreak_cutoff=50):
    import epidemic_sim
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)

    prob_outbreak = lambda outbreak_sizes: 1 - outbreak_sizes[:outbreak_cutoff-1].sum()

    if params == None:
        params  = {'verbose':False, 'edge_edit_rate':[0.05, 0.04], 'dont_cutoff_leafs':True, 'node_edit_rate':[0.05, 0.04], 'node_growth_rate':[0., 0], 'locality_bias_correction':[0.], 'enforce_connected':True, 'accept_chance_edges':1.0}
    params['preserve_degree'] = True 
    params['epsilon']         = sum(params['edge_edit_rate'])+sum(params['node_edit_rate'])
    print('Params:')
    print(params)

    stepping = 2
    max_time = 79

    if dynamics_params == None:
        #dynamics_params = {'minReplications':500, 'tie_stability':1.0, 'latent_interval':3, 'infectious_interval':300, 'max_time':250}
        #flu (per day)
        #dynamics_params = {'minReplications':1000, 'tie_stability':1.0, 'latent_interval':2, 'infectious_interval':9, 'max_time':max_time}
        dynamics_params = {'minReplications':500, 'tie_stability':1.0, 'latent_interval':3, 'infectious_interval':2, 'max_time':max_time}
    print('dynamics_params:')
    print(dynamics_params)

    if alternative_algs==None:
        alternative_algs = []
        alternative_algs += [('MUSKETEER', algorithms.generate_graph)]
        alternative_algs += [('Edge swapping', alternatives.random_noise_replicate)]
        alternative_algs += [('Expected degree model', alternatives.expected_degree_replicate)]
        alternative_algs += [('Scale-free model', alternatives.scalefree_replicate)]
        alternative_algs += [('Kronecker model', alternatives.kronecker_replicate)]
    else:
        alternative_algs = [('MUSKETEER', algorithms.generate_graph)] + alternative_algs
    vals_of_graph = []
    vals_of_replicas = {}
    for alg_name,alg_func in alternative_algs:
        vals_of_replicas[alg_name] = []

    print()
    print(G.name)

    taus = np.arange(0,1.01, 0.05)        
    for tau in taus:
        dynamics_params['tau'] = tau
        dummy, outbreak_sizes = epidemic_sim.extentSIR_long(G=G, params=dynamics_params)

        vals_of_graph.append(prob_outbreak(outbreak_sizes))

    for alg_i, (alg_name, alt_alg) in enumerate(alternative_algs):
        print()
        print(alg_name)
        for tau in taus:
            sys.stdout.write('%.2f'%tau)
            dynamics_params['tau'] = tau
            probs = []
            for replica_idx in range(num_replicas):
                replica     = alt_alg(G, params=params)
                sys.stdout.write('.')
                extent, outbreak_sizes = epidemic_sim.extentSIR_long(G=replica, params=dynamics_params)
                probs.append(prob_outbreak(outbreak_sizes))
                sys.stdout.write(' ')
                sys.stdout.flush()
            vals_of_replicas[alg_name].append(np.average(probs))
            print()

    pylab.figure()
    pylab.hold(True)

    pylab.plot(taus, vals_of_graph, '-^', color='b', linewidth=2., label='original network')

    #curvetypes = ['r-', 'g-', 'b-', 'k-', 'c-.', 'y-', 'm-D', 'b-.', 'r-.', 'c-']
    curvetypes = ['g-o', 'rd-', 'ks-', 'cp-', 'yH-', 'mD-', 'b3-', 'r4-', 'c1-']
    for alg_i, (alg_name, alt_alg) in enumerate(alternative_algs):
        print()
        print(alg_name)
        pylab.plot(taus, vals_of_replicas[alg_name], curvetypes[alg_i], linewidth=1., label=alg_name+'')

    pylab.xlabel('Transmissibility (per contact-day)', fontsize='20')
    pylab.ylabel('Probability of outbreak of size %d or more'%outbreak_cutoff, fontsize='20')#, x=0.1)
    #pylab.title(G.name)
    pylab.legend(loc='best')
    pylab.xlim(-0.01,1.01)
    pylab.ylim(ymin=-0.01)

    if figpath == None:
        figpath = 'output/alternatives_seir_outbreak_sizes_'+G.name+'_'+str(seed)+'__'+timeNow()
        figpath = clean_path(figpath)
    save_figure_helper(figpath)
    pylab.hold(False)

    data = {'vals_of_graph':vals_of_graph, 'vals_of_replicas':vals_of_replicas, 'params':params, 'dynamics_params':dynamics_params}
    graphutils.safe_pickle(path=figpath+'.pkl', data=data)


def alternatives_SEIR(G=None, seed=None, figpath=None, params=None, dynamics_params=None, alternative_algs=None, num_replicas = 150):
    import epidemic_sim
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)

    if params == None:
        params  = {'verbose':False, 'edge_edit_rate':[0.05, 0.04], 'dont_cutoff_leafs':True, 'node_edit_rate':[0.05, 0.04], 'node_growth_rate':[0., 0], 'locality_bias_correction':[0.], 'enforce_connected':True, 'accept_chance_edges':1.0}
    params['preserve_degree'] = True 
    params['epsilon']         = sum(params['edge_edit_rate'])+sum(params['node_edit_rate'])
    print('Params:')
    print(params)
    if G.name == 'potterat':
        params['matrix'] = '0.629063, 0.584033; 0.668222, 0.13873'

    stepping = 2
    max_time = 79

    if dynamics_params == None:
        #dynamics_params = {'minReplications':500, 'tau':0.1, 'tie_stability':1.0, 'latent_interval':3, 'infectious_interval':300, 'max_time':250}
        #flu (per day)
        dynamics_params = {'minReplications':1000, 'tau':0.5, 'tie_stability':1.0, 'latent_interval':2, 'infectious_interval':9, 'max_time':max_time}
        #dynamics_params = {'minReplications':100, 'tau':0.4, 'tie_stability':1.0, 'latent_interval':3, 'infectious_interval':2, 'max_time':max_time}
    print('dynamics_params:')
    print(dynamics_params)

    if alternative_algs==None:
        alternative_algs = []
        alternative_algs += [('MUSKETEER', algorithms.generate_graph)]
        alternative_algs += [('Edge swapping', alternatives.random_noise_replicate)]
        alternative_algs += [('Expected degree model', alternatives.expected_degree_replicate)]
        alternative_algs += [('Scale-free model', alternatives.scalefree_replicate)]
        alternative_algs += [('Kronecker model', alternatives.kronecker_replicate)]
    else:
        alternative_algs = [('MUSKETEER', algorithms.generate_graph)] + alternative_algs
    replicas = {'MUSKETEER':[]}
    vals_of_replicas = []    
    vals_of_replicas2 = {}
    for alg_name,alg_func in alternative_algs:
        replicas[alg_name] = []
        vals_of_replicas2[alg_name] = []
    
    mean_of_graph, outbreak_sizes = epidemic_sim.extentSIR_long(G=G, params=dynamics_params)

    print()
    print(G.name)
    pylab.figure()
    pylab.hold(True)

    xvals = list(range(dynamics_params['max_time']))
    pylab.plot(xvals[::2], mean_of_graph[::2],    '-^', color='b', linewidth=2., label='original network')
    #pylab.errorbar(xvals[::stepping], std_of_graph[::stepping], yerr=std_of_graph[::stepping], color='b', linestyle='', ecolor='b', capsize=4)[0]

    #curvetypes = ['r-', 'g-', 'b-', 'k-', 'c-.', 'y-', 'm-D', 'b-.', 'r-.', 'c-']
    curvetypes = ['g-o', 'rd-', 'ks-', 'cp-', 'yH-', 'mD-', 'b3-', 'r4-', 'c1-']
    for alg_i, (alg_name, alt_alg) in enumerate(alternative_algs):
        print()
        print(alg_name)
        for replica_idx in range(num_replicas):
            replica     = alt_alg(G, params=params)
            replicas[alg_name].append(replica)
            sys.stdout.write('.')
            extent, outbreak_sizes = epidemic_sim.extentSIR_long(G=replica, params=dynamics_params)
            sys.stdout.write(' ')
            vals_of_replicas2[alg_name].append(extent)
        vals_of_replicas2_ar = np.array(vals_of_replicas2[alg_name])
        mean_of_replicas2_ar = vals_of_replicas2_ar.mean(axis=0)
        std_of_replicas2  = vals_of_replicas2_ar.std(axis=0, ddof=1)
        print()
        pylab.plot(xvals[::2], mean_of_replicas2_ar[::2], curvetypes[alg_i], linewidth=1., label=alg_name+'')
        pylab.errorbar(xvals[::stepping], mean_of_replicas2_ar[::stepping], yerr=std_of_replicas2[::stepping], color=curvetypes[alg_i][0], linestyle='', ecolor=curvetypes[alg_i][0], capsize=2)[0]

    pylab.ylabel('New cases', fontsize='20')
    pylab.xlabel('Time (days)', fontsize='20')#, x=0.1)
    #pylab.ylabel('Normed by mean value', rotation=90, fontsize='20')
    #pylab.title(G.name)
    pylab.legend(loc='best')
    pylab.xlim(-0.01,len(mean_of_graph)+0.1)
    pylab.ylim(ymin=-0.01)

    if figpath == None:
        figpath = 'output/alternatives_seir_'+G.name+'_'+str(seed)+'__'+timeNow()
        figpath = clean_path(figpath)
    save_figure_helper(figpath)
    pylab.hold(False)

    data = {'replicas':replicas, 'vals_of_replicas':vals_of_replicas, 'vals_of_replicas2':vals_of_replicas2, 'params':params, 'dynamics_params':dynamics_params}
    graphutils.safe_pickle(path=figpath+'.pkl', data=data)

    return replicas



def bias_benchmark(G, alg, params={}, alg_params={}, seed=None):
    '''
    determine if alg has a systematic bias when running on graph G
    '''
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)
 
    if hasattr(G, 'name'):
        print('Testing bias on: ' + str(G.name))
    print('  note: error format:  NEW_VALUE - OLD_VALUE')
    if alg_params == {}:
        alg_params['verbose'] = False
        alg_params['node_edit_rate'] = [0.0]
        alg_params['edge_edit_rate'] = [0.05]
        alg_params['locality_bias_correction'] = [0.0]
        #alg_params['deferential_detachment_factor'] = 1.0
        #alg_params['edge_welfare_fraction'] = 1.0
        #alg_params['enforce_connected'] = False

    num_samples    = params.get('num_samples', 50)
    #target_metrics = params.get('target_metrics', ['num edges','clustering'])
    target_metrics = params.get('target_metrics', [metric['name'] for metric in graphutils.default_metrics if metric['optional'] == 0])

    errors_by_metric = {}
    for met_name in target_metrics:
        errors_by_metric[met_name] = []

    for it in range(num_samples):
        replica = alg(original=G, params=alg_params)

        error_dat = graphutils.compare_nets(old_G=G, new_G=replica, params={'normalize':False, 'verbose':False})[0]
        for met_name in target_metrics:
            errors_by_metric[met_name] += [error_dat.get(met_name, 0.)]
        print('Sample %d finished'%(it+1))
    print() 
    print('Analysis results')
    for met_name in errors_by_metric:
        print(met_name)
        samples = np.array(errors_by_metric[met_name])
        expected_mean = 0.
        sample_mean   = np.average(samples)
        if sum(np.abs(samples)) == 0.:
            print('metric ok.  no variation')
            print()
            continue
        try:
            t_val, p_val  = scipy.stats.ttest_1samp(a=samples, popmean=expected_mean)
        except:
            print('Couldn\'t compute: '+met_name)
            print()
            continue
        print('    Mean error: %.9f, Expected: %.9f'%(sample_mean,expected_mean))
        print('    p_val:      %.9f'%(p_val,))
        print('    %s at 5%% level'%(p_val<0.05 and 'REJECTED' or 'apparently OK'))
        print(samples)
        print()
 
    return errors_by_metric

def clean_path(fpath):
#make the path suitable for TeX
    fpath = fpath.replace('.', 'd')
    fpath = fpath.replace(',', 'c')
    fpath = fpath.replace(' ', '_')
    fpath = fpath.replace('*', '_')
    fpath = fpath.replace('(', '_')
    fpath = fpath.replace(')', '_')

    return fpath


def compare_nets_betweenness(old_G, new_G, params=None):
    '''
    Report on the differences between two networks for their betweenness
    '''
    print('Node\tBC\tNewBC\tDIFF\%\tDC\tNewDC\tDIFF\%')
    old_bet_cet = nx.betweenness_centrality(old_G)
    new_bet_cet = nx.betweenness_centrality(new_G)
    bet_errors  = []
    deg_errors  = []
    for node in old_G:  #assumes same nodes
        bet_old = old_bet_cet[node]
        bet_new = new_bet_cet[node]
        bet_error = (bet_new-bet_old)/max(bet_old,bet_new)
        deg_old = old_G.degree(node)
        deg_new = new_G.degree(node)
        deg_error = (deg_new-deg_old)/max(deg_old,deg_new)
        print('%s\t%.2f\t%.2f\t%.2f%%\t%.2f\t%.2f\t%.2f%%'%(str(node),bet_old,bet_new,100*bet_error,deg_old,deg_new,100*deg_error))
        bet_errors.append(bet_error)
        deg_errors.append(deg_error)
    print('Mean betweenness difference: %.2f%%'%np.average(bet_errors))
    print('Mean degree      difference: %.2f%%'%np.average(deg_errors))


def degree_expect_diff(old_G, new_G):
    #the expectation of d over degree histograph of old_G  
    #   minus the expectation of d over degree histograph of new_G
    #(normalized by the largest expectation)
    h1 = nx.degree_histogram(old_G)
    h2 = nx.degree_histogram(new_G)
    if len(h1) < len(h2):
        h1 += [0]*(len(h2)-len(h1))
    else:
        h2 += [0]*(len(h1)-len(h2))
    ret = 0.
    denom_1 = 0.
    denom_2 = 0.
    for i,val2 in enumerate(h2):
        ret += i*(val2 - h1[i])
        denom_1 += i*h1[i]
        denom_2 += i*val2
    return ret/max(denom_1,denom_2)


def default_graphs():
    print('Loading default graphs:')
    graphs = {}
    names = [('../benchmark/net4-1.graph',         'adjlist',  'engineering',),
             ]

    for graph_name, file_type, graph_type in names:
        G = graphutils.load_graph(path=graph_name, params={'graph_type':file_type})
        graphs[graph_name] = {'G':G, 'graph_type':graph_type}
        print(graph_name)

    graphs['southern women'] = {'G':nx.generators.davis_southern_women_graph(), 'graph_type':'social'}
    print('.')
    graphs['karate club']    = {'G':nx.generators.karate_club_graph(),          'graph_type':'social'}
    print('.')
    
    graphs['erdos-renyi200'] = {'G':nx.erdos_renyi_graph(n=200, p=0.02, seed=42),          'graph_type':'classic_model'}
    print('.')
    graphs['watts-strogatz200'] = {'G':nx.watts_strogatz_graph(n=200, k=3, p=0.05, seed=42),'graph_type':'classic_model'}
    print('.')
    graphs['barabasi-albert200'] = {'G':nx.barabasi_albert_graph(n=200, m=10, seed=42),'graph_type':'classic_model'}
    print('.')
    random.seed()

    while True:
        yield graphs

def dissimilarity_preservation_demo():
    metrics = graphutils.default_metrics[:]
    metrics = [met for met in metrics if met['name'] not in ['avg flow closeness', 'powerlaw exp', 'density', 'num comps', 'harmonic mean path', 'mean ecc', 'average degree']]
    metrics = [met for met in metrics if met['name'] not in ['degree assortativity']]  #varies too much
    metrics = [met for met in metrics if met['optional'] == 0]

    ER = nx.erdos_renyi_graph(n=100, p=0.08)
    ER.name = 'ER'
    scaling_data, norms = self_dissimilarity_analysis(G=ER, metrics=metrics, markersize=20)

    params  = {'verbose':False, 'edge_edit_rate':[0.03], 'node_edit_rate':[], 'node_growth_rate':[0], 
            'enforce_connected':True, 'accept_chance_edges':1.0, 'locality_bias_correction':[-0.65, 0, 0, 0, 0]}

    for rep_num in range(10):
        replica = algorithms.generate_graph(original=ER, params=params)

        self_dissimilarity_analysis(G=replica, metrics=metrics, norms = norms, figure=1, markersize=1)
    
    pylab.savefig('output/self_asymmetry_compare__'+timeNow()+'.eps')


def editing_demo(seed=15 ):
    original = graphutils.load_graph('data-cyber-small/gr2.gml')
    #maybe: convert_to_integers to simplify display?
    base_params      = {'edge_edit_rate':[], 'node_edit_rate':[], 'node_growth_rate':[]}
    
    npr.seed(seed)
    random.seed(seed)
    pos = nx.spring_layout(original) 
    
    
    editing_demo_draw(G=original, new_G=original, pos=pos, seed=1)

    npr.seed(15)
    random.seed(15)
    params_blue_edge = base_params.copy()
    params_blue_edge['edge_edit_rate'] = [0.05]
    blue_edge_replica = algorithms.generate_graph(original=original, params=params_blue_edge)
    editing_demo_draw(G=original, new_G=blue_edge_replica, pos=pos, seed=10)
    
    npr.seed(62)
    random.seed(62)
    params_red_edge = base_params.copy()
    params_red_edge['edge_edit_rate'] = [0, 0.2]
    red_edge_replica = algorithms.generate_graph(original=original, params=params_red_edge)
    editing_demo_draw(G=original, new_G=red_edge_replica, pos=pos, seed=18    )
    
    npr.seed(    2    )
    random.seed(  2   )
    params_blue_node = base_params.copy()
    params_blue_node['node_edit_rate'] = [0.05]
    blue_node_replica = algorithms.generate_graph(original=original, params=params_blue_node)
    editing_demo_draw(G=original, new_G=blue_node_replica, pos=pos, seed=22)
    

    npr.seed(    36  )  #7,18 ,24
    random.seed(  36 )
    params_red_node = base_params.copy()
    params_red_node['node_edit_rate'] = [0, 0.2]
    red_node_replica = algorithms.generate_graph(original=original, params=params_red_node)
    editing_demo_draw(G=original, new_G=red_node_replica, pos=pos, seed=9)

def editing_demo_draw(G, new_G, seed=1, **kwargs):
    #prepare a graph
    npr.seed(seed)
    random.seed(seed)

    delta = graphutils.graph_graph_delta(G, new_G)
    new_nodes = delta['new_nodes']
    del_nodes = delta['del_nodes']
    new_edges = delta['new_edges']
    del_edges = delta['del_edges']

    if 'pos' not in kwargs:
        merged_pos = nx.spring_layout(G)
    else:
        merged_pos = kwargs['pos'].copy()

    merged_G = G.copy()
    tmp_pos  = merged_pos.copy()
    for node in new_nodes:
        merged_G.add_node(node)
        tmp_pos[node] = npr.rand(2)
    for edge in new_edges:
        merged_G.add_edge(*edge)
    tmp_pos = nx.spring_layout(merged_G, pos=tmp_pos)
    for node in new_nodes:
        merged_pos[node]=tmp_pos[node]

    regular_nodes = [node for node in merged_G         if ((node not in new_nodes) and (node not in del_nodes))]
    regular_edges = [edge for edge in merged_G.edges() if (G.has_edge(*edge) and new_G.has_edge(*edge))]
    labels = {}
    for node in merged_G:
        labels[node] = ''


    base_node_color = 'black'
    new_node_color  = 'brown'
    del_node_color  = 'yellow' #'blue'
    base_edge_color = 'black'
    new_edge_color  = new_node_color
    del_edge_color  = 'black'
   
    base_edge_style = 'solid'
    new_edge_style  = 'solid'
    del_edge_style  = 'dashed'
    new_edge_width  = 5
    del_edge_width  = 2

    #pylab.close()
    #pylab.figure(figsize=(2.1, 2))
    pylab.figure()
    pylab.hold(True)
    #nx.draw(G)
    #nx.draw(G, pos=pos, labels=labels, width=widths,  style=styles)
    nx.draw_networkx_nodes(merged_G, pos=merged_pos, nodelist=regular_nodes, node_color=base_node_color, node_size=500, alpha=1.0, with_labels=False, node_shape='s')
    nx.draw_networkx_nodes(merged_G, pos=merged_pos, nodelist=new_nodes, node_color=new_node_color, node_size=500, alpha=1.0, with_labels=False, node_shape='s')
    nx.draw_networkx_nodes(merged_G, pos=merged_pos, nodelist=del_nodes, node_color=del_node_color, node_size=100, alpha=1.0,  with_labels=False, node_shape='s')
    nx.draw_networkx_edges(merged_G, pos=merged_pos, edgelist=regular_edges, style=base_edge_style, edge_color=base_edge_color, alpha=1.0)  #width=dict
    nx.draw_networkx_edges(merged_G, pos=merged_pos, edgelist=new_edges, style=new_edge_style, edge_color=new_edge_color, alpha=1.0, width=new_edge_width)
    nx.draw_networkx_edges(merged_G, pos=merged_pos, edgelist=del_edges, style=del_edge_style, edge_color=del_edge_color, alpha=1.0, width=del_edge_width)
    nx.draw_networkx_labels(G, pos=merged_pos, labels=labels)
    pylab.grid(b=False)
    #pylab.axis('on')
    pylab.axis('off')
    #pylab.axes(frameon=False)

   
    #pylab.savefig('output/graph_delta_'+timeNow() + '.pdf')
    pylab.savefig('output/graph_delta_'+timeNow() + '.eps')
    pylab.hold(False)
    time.sleep(2)

    #pylab.show()
   
    npr.seed()
    random.seed()
    
    return delta 

def ensemble_validation_benchmark(alg, params={}, graphs={}, alg_params={}):
    '''
    compare a population the replicas against a population of originals
    '''
    raise NotImplemented


    if graphs == {}:
        graphs = next(default_graphs())
    if alg_params == {}:
        alg_params['verbose'] = False
        alg_params['error_rate'] = [0.05/(1.+i) for i in range(100)]
    if params == {}:
        params['num_replicas'] = 10

    replica_ensembles = {}
    for graph_name in graphs:
        G          = graphs[graph_name]['G']
        graph_type = graphs[graph_name]['graph_type']

        replica_ensembles[graph_name] = []

        for r in range(params['num_replicas']):
            replica = alg(original=G, params=alg_params)
            replica_ensembles[graph_name].append(replica)

    #TODO: plot the originals and their replicas over several metrics


def movie_of_replicas():
    base_graph = graphutils.load_graph('data-social/911_unwted.gml')
    base_graph.name = '911'
    #base_graph = graphutils.load_graph(path='data-epi/potterat.gml')
    #base_graph.name = 'potterat'
    #assert base_graph.number_of_nodes() == 250
    #assert base_graph.number_of_edges() == 266
    num_replicas = 500
    #base_graph = graphutils.load_graph(path='data-social/911_unwted.gml')
    #base_graph.name = '911_krebs'
    params_potterat  = {'verbose':False, 'edge_edit_rate':[0.08, 0.07], 'node_edit_rate':[0.08, 0.07], 'node_growth_rate':[0], 'enforce_connected':True, 'accept_chance_edges':1.0, 'new_edge_horizon':20, 'num_deletion_trials':20, 'locality_bias_correction':[0, 0, 0, 0, 0]}

    params_911 = {'edge_edit_rate':[0.08, 0.07], 'node_edit_rate':[0.08, 0.07], 'verbose':False, 'new_edge_horizon':20, 'num_insertion_trials':30, 'enforce_connected':True, 'locality_bias_correction':[1.0, 1.0, 0, 0, 0], 'edit_method':'alternating', 'deferential_detachment_factor':1.,  'edge_welfare_fraction':1.0}

    try:
        os.mkdir('output/movie')
    except:
        pass
    gpath     = 'output/movie/'+base_graph.name+'_'+timeNow()+'.dot'
    gpath_fig = gpath[:-3]+'png'
    graphutils.write_graph(G=base_graph, path=gpath)
    visualizer_cmdl = 'sfdp -Nlabel="" -Nfontsize=0 -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Nfillcolor=/ylgnbu9/9 \
                            -Nstyle=filled -Nperipheries=0 -Gratio=0.75 -Gdpi=300 -Tpng %s > %s &'%(gpath,gpath_fig)
    #http://www.graphviz.org/doc/info/colors.html
    print('Writing graph image: %s ..'%gpath_fig)
    retCode = os.system(visualizer_cmdl)
    
    for replica_num in range(num_replicas):
        #gpath_rep = gpath[:-4]+'_replica%d'%(10000+replica_num)+'.dot'
        gpath_rep = gpath[:-4]+'_replica%d'%(replica_num)+'.dot'
        gpath_rep_fig = gpath_rep[:-4]+'.png'
        replica = algorithms.generate_graph(original=base_graph, params=params_911)
        graphutils.write_graph(G=replica, path=gpath_rep)
        visualizer_cmdl = 'sfdp -Nlabel="" -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Nfillcolor=/ylgnbu9/9 \
                                -Nstyle=filled -Nperipheries=0 -Gratio=0.75 -Gdpi=300 -Tpng %s > %s &'%(gpath_rep,gpath_rep_fig)
        #http://www.graphviz.org/doc/info/colors.html
        print('Writing replica image: %s ..'%gpath_rep_fig)
        retCode = os.system(visualizer_cmdl)


def paper_benchmarks_and_figures():
#    execute all the benchmarks used in the paper
    paper_findings_0()
    paper_findings_1()
    paper_findings_2()
    paper_findings_emergent()
    paper_findings_one_metric()
    paper_findings_snapshots()
    paper_alternatives_benchmarks()
    paper_findings_memoriless()
    paper_findings_epidemic()

    paper_illustration_appendix()

def paper_alternatives_benchmarks():
    base_graph = graphutils.load_graph(path='data-epi/potterat.gml')
    base_graph.name = 'potterat'
    num_replicas = 150

    metrics = graphutils.default_metrics[:]
    metrics  = [met for met in metrics if met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree connectivity', 'degree assortativity', ]]
    params_MUSKETEER  = {'verbose':False, 'edge_edit_rate':[0, 0, 0,  0, 0, 0.08], 'node_edit_rate':[0, 0, 0,  0, 0, 0.08], 'node_growth_rate':[0., 0], 'enforce_connected':True, 'accept_chance_edges':1.0}
    params_rand_noise = {'epsilon':sum(params_MUSKETEER['edge_edit_rate'])+sum(params_MUSKETEER['node_edit_rate']), 'preserve_degree':False}
    params_edge_swap = params_rand_noise.copy()
    params_edge_swap['preserve_degree'] = True
    params_edge_swap['preserve_connected'] = True

    #print 'Musketeer:'
    testers.replica_vs_original(G=base_graph, num_replicas=num_replicas, seed=None, params=params_MUSKETEER)
    #print

    '''
    print 'ER:'
    testers.replica_vs_original(G=base_graph, num_replicas=num_replicas, seed=None, params={}, metrics=metrics, 
            generator_func=alternatives.er_replicate, title_infix='er')
    
    print
    print 'scalefree:'
    testers.replica_vs_original(G=base_graph, num_replicas=num_replicas, seed=None, params={}, metrics=metrics,
            generator_func=alternatives.scalefree_replicate, title_infix='scale_free')
    print
    print 'smallworld:'
    testers.replica_vs_original(G=base_graph, num_replicas=num_replicas, seed=None, params={}, metrics=metrics,                 
            generator_func=alternatives.watts_strogatz_replicate, title_infix='small_world')
    print 
    print 'Expected degree:'
    testers.replica_vs_original(G=base_graph, num_replicas=num_replicas, seed=None, params={}, metrics=metrics,
            generator_func=alternatives.expected_degree_replicate, title_infix='expected_deg')
    
    print
    print 'random noising:'
    testers.replica_vs_original(G=base_graph, num_replicas=num_replicas, seed=None, params=params_rand_noise, metrics=metrics,
            generator_func=alternatives.random_noise_replicate, title_infix='rand_edit')
    print
    print 'edge swap:'
    testers.replica_vs_original(G=base_graph, num_replicas=num_replicas, seed=None, params=params_edge_swap, metrics=metrics, 
            generator_func=alternatives.random_noise_replicate, title_infix='edge_swap')
    '''
    print()
    print('kronecker:')
    #FIXME kron_fitted_matrix = alternatives.kronecker_replicate(original=base_graph, params={'just_do_fitting':True})
    print('Kronecker fitted matrix:')
    #kron_fitted_matrix = '0.632785, 0.583371; 0.666826, 0.14198'
    kron_fitted_matrix = '0.629063, 0.584033; 0.668222, 0.13873'
    print(kron_fitted_matrix)
    params_kronecker = {'matrix':kron_fitted_matrix}
    testers.replica_vs_original(G=base_graph, num_replicas=num_replicas, seed=None, params=params_kronecker, metrics=metrics,
            generator_func=alternatives.kronecker_replicate, title_infix='kronecker')
    
def paper_findings_0(seed=10):
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)

    base_graph = graphutils.load_graph(path='data-samples/mesh33.edges')
    print('Musketeer:')
    params_MUSKETEER  = {'edge_edit_rate':[0.01]}
    
    ######################
    #     1st figure
    ######################
    #http://www.graphviz.org/doc/info/colors.html
    gpath     = 'output/'+base_graph.name+'_'+timeNow()+'.dot'
    gpath_fig = gpath[:-3]+'eps'
    graphutils.write_graph(G=base_graph, path=gpath)
    #visualizer_cmdl = 'sfdp -Nlabel="" -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Nfillcolor=/ylgnbu9/9 -Nstyle=filled -Nperipheries=0 -Teps %s > %s &'%(gpath,gpath_fig)
    visualizer_cmdl = 'sfdp  -Nlabel="" -Nwidth=0.06 -Nfixedsize=true -Nheight=0.06 -Teps %s > %s &'%(gpath,gpath_fig)
    print('Writing graph image: %s ..'%gpath_fig)
    retCode = os.system(visualizer_cmdl)
    
    gpath_rep = gpath[:-4]+'_replica_shallow.dot'
    gpath_rep_fig = gpath_rep[:-4]+'.eps'
    params2 = params_MUSKETEER.copy()
    params2['verbose']=True
    replica = algorithms.generate_graph(original=base_graph, params=params2)
    graphutils.write_graph(G=replica, path=gpath_rep)
    visualizer_cmdl = 'sfdp  -Nlabel="" -Nwidth=0.06 -Nfixedsize=true -Nheight=0.06 -Teps %s > %s &'%(gpath_rep,gpath_rep_fig)
    print('Writing graph image: %s ..'%gpath_rep_fig)
    retCode = os.system(visualizer_cmdl)

    ######################
    #     2nd figure
    ######################
    gpath_rep = gpath[:-4]+'_replica_deep.dot'
    gpath_rep_fig = gpath_rep[:-4]+'.eps'
    params3 = params_MUSKETEER.copy()
    params3['edge_edit_rate'] = [0, 0, 0, 0, 0.01]
    replica = algorithms.generate_graph(original=base_graph, params=params3)
    graphutils.write_graph(G=replica, path=gpath_rep)
    visualizer_cmdl = 'sfdp  -Nlabel="" -Nwidth=0.06 -Nfixedsize=true -Nheight=0.06 -Teps %s > %s &'%(gpath_rep,gpath_rep_fig)
    print('Writing graph image: %s ..'%gpath_rep_fig)
    retCode = os.system(visualizer_cmdl)



def paper_findings_1(seed=8, intermediates=False, num_replicas=150):
    #3 is good
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)

    base_graph = graphutils.load_graph(path='data-epi/potterat.gml')
    base_graph.name = 'potterat'
    assert base_graph.number_of_nodes() == 250
    assert base_graph.number_of_edges() == 266
    #base_graph = graphutils.load_graph(path='data-social/911_unwted.gml')
    #base_graph.name = '911_krebs'
    
    print('Musketeer:')
    params_MUSKETEER  = {'verbose':False, 'edge_edit_rate':[0, 0, 0, 0, 0.08], 'node_edit_rate':[0, 0, 0, 0, 0.08], 'node_growth_rate':[0], 'dont_cutoff_leafs':True,
                        'enforce_connected':True, 'accept_chance_edges':0.8, 'new_edge_horizon':20, 'num_deletion_trials':20, 'locality_bias_correction':[0,0]}
    metrics = graphutils.default_metrics[:]
    metrics  = [met for met in metrics if met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree connectivity', 'degree assortativity', ]]
    #metrics  = filter(lambda met: met['name'] not in ['avg harmonic shortest path'], metrics)

    
    ######################
    #     1st figure
    ######################
    #http://www.graphviz.org/doc/info/colors.html
    gpath     = 'output/'+base_graph.name+'_'+timeNow()+'.dot'
    gpath_fig = gpath[:-3]+'eps'
    graphutils.write_graph(G=base_graph, path=gpath)
    #visualizer_cmdl = 'sfdp -Nlabel="" -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Nfillcolor=/ylgnbu9/9 -Nstyle=filled -Nperipheries=0 -Teps %s > %s &'%(gpath,gpath_fig)
    visualizer_cmdl = 'sfdp  -Nlabel="" -Nwidth=0.06 -Nfixedsize=true -Nheight=0.06 -Teps %s > %s &'%(gpath,gpath_fig)
    print('Writing graph image: %s ..'%gpath_fig)
    retCode = os.system(visualizer_cmdl)
    
    gpath_rep = gpath[:-4]+'_replica.dot'
    gpath_rep_fig = gpath_rep[:-4]+'.eps'
    params2 = params_MUSKETEER.copy()
    params2['verbose']=True
    replica = algorithms.generate_graph(original=base_graph, params=params2)
    graphutils.write_graph(G=replica, path=gpath_rep)
    visualizer_cmdl = 'sfdp  -Nlabel="" -Nwidth=0.06 -Nfixedsize=true -Nheight=0.06 -Teps %s > %s &'%(gpath_rep,gpath_rep_fig)
    print('Writing graph image: %s ..'%gpath_rep_fig)
    retCode = os.system(visualizer_cmdl)

    ######################
    #     2nd figure
    ######################
    #testers.replica_vs_original(G=base_graph, num_replicas=num_replicas, seed=1, params=params_MUSKETEER, metrics=metrics, title_infix='musketeer', intermediates=intermediates)

    replica_vs_original_degree_distribution(G=base_graph, num_replicas=num_replicas, seed=1, params=params_MUSKETEER)
    


def paper_findings_2():
    num_replicas = 150
    params_MUSKETEER  = {'verbose':False, 'edge_edit_rate':[0, 0, 0, 0.08], 'node_edit_rate':[0, 0, 0, 0.08], 'node_growth_rate':[0], 
            'enforce_connected':True, 'accept_chance_edges':1.0, 'num_pairs_to_sample':100, 'num_trial_particles':50, 'num_insertion_trials':50}
    metrics = graphutils.default_metrics[:]
    metrics  = [met for met in metrics if met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree assortativity', 'degree connectivity', 'powerlaw exp']]

    ######################
    #     3rd figure
    ######################
    print('Barabasi-Albert')
    params_BA = params_MUSKETEER.copy()
    params_BA['locality_bias_correction'] = []
    metrics += [met for met in graphutils.default_metrics if met['name'] == 'powerlaw exp']
    metrics[-1]['optional']=0
    BA_graph = nx.barabasi_albert_graph(n=300, m=10, seed=42)
    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_BA, metrics=metrics)

    print()
    print('Erdos-Renyi')
    params_ER = params_MUSKETEER.copy()
    params_ER['locality_bias_correction'] = [] 
    ER_graph = nx.erdos_renyi_graph(n=300, p=0.04, seed=42)
    metrics_ER = [met for met in metrics if met['name'] not in ['powerlaw exp']]
    testers.replica_vs_original(G=ER_graph, num_replicas=num_replicas, seed=1, params=params_ER, metrics=metrics_ER)

    ##WS_graph = nx.watts_strogatz_graph(n=300, k=4, p=0.001)
    ##metrics_WS = filter(lambda met: met['name'] not in ['powerlaw exp', 'degree assortativity'], metrics)
    ##testers.replica_vs_original(G=WS_graph, num_replicas=num_replicas, seed=1, params=params_MUSKETEER, metrics=metrics_WS)


def paper_findings_emergent():
    base_graph = graphutils.load_graph(path='data-epi/potterat.gml')
    base_graph.name = 'potterat'
    assert base_graph.number_of_nodes() == 250
    assert base_graph.number_of_edges() == 266
    num_replicas = 150
    #base_graph = graphutils.load_graph(path='data-social/911_unwted.gml')
    #base_graph.name = '911_krebs'
    
    print('Musketeer:')
    params_MUSKETEER  = {'verbose':False, 'edge_edit_rate':[0, 0, 0, 0, 0, 0.08], 'node_edit_rate':[0, 0, 0, 0, 0, 0.08], 'node_growth_rate':[0], 'dont_cutoff_leafs':True, 'enforce_connected':True, 'accept_chance_edges':0.8, 'new_edge_horizon':20, 'num_deletion_trials':20, 'locality_bias_correction':[0, 0, 0, 0, 0],
    'skip_param_sanity_check':True}

    
    ######################
    #     5a figure
    ######################
    replicas = alternatives_SEIR(G=base_graph, num_replicas=num_replicas, seed=1, params=params_MUSKETEER)

    ######################
    #     5b figure
    ######################
    #alternatives_outbreaks(G=base_graph, num_replicas=num_replicas, seed=1, params=params_MUSKETEER)
    

def paper_findings_epidemic(seed=8):
    #3 is good
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)

    num_replicas = 20
    #num_replicas = 150
    
    params_default  = {'verbose':False, 'edge_edit_rate':[0.08, 0.07], 'node_edit_rate':[0.08, 0.07], 'node_growth_rate':[0], 
            'dont_cutoff_leafs':False, #NOTE: this is important
            'new_edge_horizon':10, 'num_deletion_trials':20, 'locality_bias_correction':[0,], 'edit_method':'sequential',
            #'accept_chance_edges':1.0, 
            #'edit_method':'alternating',#NOTE: this is important
            #'enforce_connected':True, 
            'enforce_connected':False, 
            #'component_is_edited': [True] + 1000*[False],
            #'component_is_edited': [False, False, False, False, True] + 1000*[True],
            #'deferential_detachment_factor':0.0, #1.0,
            }
    #params_default['algorithm'] = algorithms.musketeer_on_subgraphs

    metrics_default = graphutils.default_metrics[:]
    metrics_default  = [met for met in metrics_default if met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree connectivity', 'degree assortativity',  'average shortest path', 'mean ecc', 'powerlaw exp', ]]
    #metrics  = filter(lambda met: met['name'] not in ['avg harmonic shortest path'], metrics)

    problems = [('adol.edgelist',params_default,metrics_default),  
                #('salganik.edgelist',params_default,metrics_default),  

                 #('clus.edgelist',params_default,metrics_default),  
                 #('haar.edgelist',params_default,metrics_default),  
                 #('urb1.edgelist',params_default,metrics_default),  
                 #('anti.edgelist',params_default,metrics_default),  
                 #('flag.edgelist',params_default,metrics_default),  
                 #('proj.edgelist',params_default,metrics_default),  
                 #('urb2.edgelist',params_default,metrics_default),
               ]

    for base_graph_name,params,metrics in problems:
        base_graph = graphutils.load_graph(path='data-epi/'+base_graph_name)
        base_graph.name = base_graph_name

        gpath     = 'output/'+base_graph.name+'_'+timeNow()+'.dot'
        gpath_fig = gpath[:-3]+'eps'
        graphutils.write_graph(G=base_graph, path=gpath)
        visualizer_cmdl = 'sfdp  -Nlabel="" -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Teps %s > %s &'%(gpath,gpath_fig)
        print('Writing graph image: %s ..'%gpath_fig)
        retCode = os.system(visualizer_cmdl)
        
        #gpath_rep = gpath[:-4]+'_replica.dot'
        #gpath_rep_fig = gpath_rep[:-4]+'.eps'
        #params2 = params.copy()
        #params2['verbose']=True
        #replica = algorithms.generate_graph(original=base_graph, params=params2)
        #graphutils.write_graph(G=replica, path=gpath_rep)
        #visualizer_cmdl = 'sfdp  -Nlabel="" -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Teps %s > %s &'%(gpath_rep,gpath_rep_fig)
        #print 'Writing graph image: %s ..'%gpath_rep_fig
        #retCode = os.system(visualizer_cmdl)

        testers.replica_vs_original(G=base_graph, num_replicas=num_replicas, seed=1, params=params, metrics=metrics, title_infix='musketeer')
        #replica_vs_original_degree_distribution(G=base_graph, num_replicas=num_replicas, seed=1, params=params_MUSKETEER)
    

def paper_findings_extended(seed=8):
#figures for the appendix
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)

    num_replicas = 30 #FIXME
    #um_replicas = 150
    
    params_default  = {'verbose':False, 
    'edge_edit_rate':[0.04, 0.03], 'node_edit_rate':[0.04, 0.03], 'node_growth_rate':[0],  
    'dont_cutoff_leafs':False, #NOTE: this is important
            'enforce_connected':True, 'accept_chance_edges':1.0, 'new_edge_horizon':10, 'num_deletion_trials':20, 'locality_bias_correction':[0,],
            'num_v_cycles':10,
            }
    params_default['algorithm'] = (algorithms.musketeer_on_subgraphs, algorithms.musketeer_iterated_cycle,)
    params_for_big = params_default.copy()
    params_for_big['new_edge_horizon'] = 3
    params_for_big['num_v_cycles'] = 1
    params_for_big['algorithm'] = algorithms.musketeer_on_subgraphs
    metrics_for_big = [k for k in graphutils.default_metrics if k['name'] in ['num nodes', 'num edges', 'clustering', 'total deg*deg', 'harmonic mean path']]

    metrics_default = graphutils.default_metrics[:]
    metrics_default  = [met for met in metrics_default if met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree connectivity', 'degree assortativity',  'average shortest path', 'mean ecc', 'powerlaw exp', ]]
    #metrics  = filter(lambda met: met['name'] not in ['avg harmonic shortest path'], metrics)

    problems = [#{'graph_path':'data/mesh33.gml',   'graph_params':params_default, 'graph_metrics':metrics_default,},
                #{'graph_path':'data-appendix-new/celegansneural_undirected_graph.elist', },
                #{'graph_path':'data-appendix-new/911_unwted.gml',     },
                #huge {'graph_path':'data-social/actor.adjlist',     },
                #{'graph_path':'data-appendix-new/arenas_email.edges',     },
                #{'graph_path':'data-appendix-new/dolphins.gml',     },
                #{'graph_path':'data-appendix-new/newman06_netscience.gml',},
                #{'graph_path':'data-appendix-new/watts_strogatz98_power.elist', },
                #{'graph_path':'data-appendix-new/adol.edgelist'},
                #{'graph_path':'data-appendix-new/salganik.edgelist'},
                #{'graph_path':'data-appendix-new/clus.edgelist'},  
                #{'graph_path':'data-appendix-new/haar.edgelist'},  
                #{'graph_path':'data-appendix-new/urb1.edgelist'},  
                #{'graph_path':'data-appendix-new/anti.edgelist'},  
                #{'graph_path':'data-appendix-new/flag.edgelist'},  
                #{'graph_path':'data-appendix-new/proj.edgelist'},  
                #{'graph_path':'data-appendix-new/urb2.edgelist'},
                ##big ones
                #{'graph_path':'data-appendix-new/HBCI-network.elist',   'graph_params':params_for_big, 'graph_metrics':metrics_for_big,},
                #{'graph_path':'data-appendix-new/dictionary28_cc1.elist',   'graph_params':params_for_big, 'graph_metrics':metrics_for_big,},
                #{'graph_path':'data-appendix-new/Lederberg_cc1.elist',   'graph_params':params_for_big, 'graph_metrics':metrics_for_big,},
                #{'graph_path':'data-appendix-new/p2p-Gnutella09.elist', 'graph_params':params_for_big, 'graph_metrics':metrics_for_big,},
                {'graph_path':'data-engineering/USairport_2010.elist', 'graph_params':params_for_big, 'graph_metrics':metrics_for_big,},
                {'graph_path':'data-social/Wiki-Vote.elist', 'graph_params':params_for_big, 'graph_metrics':metrics_for_big,},
                {'graph_path':'data-engineering/as-caida20071105.welist', 'graph_params':params_for_big, 'graph_metrics':metrics_for_big,},
                {'graph_path':'data-social/OClinks_w_chars.welist', 'graph_params':params_for_big, 'graph_metrics':metrics_for_big,},
                {'graph_path':'data-social/PGPgiantcompo.net', 'graph_params':params_for_big, 'graph_metrics':metrics_for_big,},
                {'graph_path':'data-biological/salmonela_barabasi__TY.elist',},
                {'graph_path':'data-biological/baydry.net',},
                {'graph_path':'data-biological/bo_yeast_barabasi.elist',},
                {'graph_path':'data-social/norwegian_boarddirectors_net1m_2011-08-01.elist',},
               ]

    for problem_info in problems:
        graph_path      = problem_info['graph_path']
        graph_params    = problem_info.get('graph_params', params_default)
        graph_metrics   = problem_info.get('graph_metrics', metrics_default)
        graph_num_replicas = problem_info.get('num_replicas', num_replicas)
        #graph_num_replicas = 2
        #FIXME print 'Starting:'
        #FIXME print graph_path
        #FIXME base_graph      = graphutils.load_graph(path=graph_path)
        base_graph      = graphutils.load_graph(path=graph_path, params={'verbose':False})
        #print 'Type: ' + str(type(base_graph))
        print('Size nn: %d, ne: %d'%(nx.number_of_nodes(base_graph), nx.number_of_edges(base_graph)))
        base_graph = nx.Graph(base_graph)
        #FIXME 
        print(graph_path, nx.number_of_nodes(base_graph), nx.number_of_edges(base_graph))
        continue

        base_graph.name = graph_path.split('/')[1]

        gpath     = 'output/'+base_graph.name+'_'+timeNow()+'.dot'
        gpath_fig = gpath[:-3]+'eps'
        graphutils.write_graph(G=base_graph, path=gpath)
        #try:
        #    visualizer_cmdl = 'sfdp  -Nlabel="" -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Teps %s > %s &'%(gpath,gpath_fig)
        #    print 'Writing graph image: %s ..'%gpath_fig
        #    retCode = os.system(visualizer_cmdl)
        #except Exception,inst:
        #    print 'Warning: could not plot '+graph_path + ': '+str(inst)   
        #    exc_traceback = sys.exc_info()[2]
        #    print str(inst) + "\n" + str(traceback.format_tb(exc_traceback)).replace('\\n', '\n')


        #gpath_rep = gpath[:-4]+'_replica.dot'
        #gpath_rep_fig = gpath_rep[:-4]+'.eps'
        #params2 = params.copy()
        #params2['verbose']=True
        #replica = algorithms.generate_graph(original=base_graph, params=params2)
        #graphutils.write_graph(G=replica, path=gpath_rep)
        #visualizer_cmdl = 'sfdp  -Nlabel="" -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Teps %s > %s &'%(gpath_rep,gpath_rep_fig)
        #print 'Writing graph image: %s ..'%gpath_rep_fig
        #retCode = os.system(visualizer_cmdl)
    
        try:
            #runs into trouble ...
            #for param_set in testers.param_set_generator(fixed_set=graph_params):
            #    testers.replica_vs_original(G=base_graph, num_replicas=graph_num_replicas, seed=seed, params=param_set, metrics=graph_metrics, title_infix='', n_jobs=-1)
            for param_set in testers.param_set_generator():
                testers.replica_vs_original(G=base_graph, num_replicas=graph_num_replicas, seed=seed, params=param_set, metrics=graph_metrics, title_infix='', n_jobs=-1)
            #replica_vs_original_degree_distribution(G=base_graph, num_replicas=num_replicas, seed=1, params=params_MUSKETEER)
            #new_G = algorithms.generate_graph(original=base_graph, params=graph_params)
            #graphutils.compare_nets(old_G=base_graph, new_G=new_G, params=graph_params, metrics=graph_metrics)
        except Exception as inst:
            print('Warning: could not generate or evaluate '+graph_path + ': '+str(inst))
            exc_traceback = sys.exc_info()[2]
            print(str(inst) + "\n" + str(traceback.format_tb(exc_traceback)).replace('\\n', '\n'))

def paper_findings_memoriless(seed=1):
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)

    base_graph = graphutils.load_graph(path='data-epi/potterat.gml')
    base_graph.name = 'potterat'
    assert base_graph.number_of_nodes() == 250
    assert base_graph.number_of_edges() == 266
    num_replicas = 150
    #base_graph = graphutils.load_graph(path='data-social/911_unwted.gml')
    #base_graph.name = '911_krebs'
    
    print('Musketeer:')
    params_MUSKETEER  = {'verbose':False, 'edge_edit_rate':[0.08, 0.07, 0.001, 0.001, 0.001, ], 'node_edit_rate':[0.08, 0.07, 0.001, 0.001, 0.001, ], 'node_growth_rate':[0], 'dont_cutoff_leafs':True,
            'enforce_connected':True, 'accept_chance_edges':1.0, 'new_edge_horizon':20, 'num_deletion_trials':20, 
             'locality_bias_correction':[-0.0, -0.2, 0, 0, 0],
            'memoriless_interpolation':True} #NOTE the complication
    metrics = graphutils.default_metrics[:]
    metrics  = [met for met in metrics if met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree connectivity', 'degree assortativity', ]]
    #metrics  = filter(lambda met: met['name'] not in ['avg harmonic shortest path'], metrics)

    testers.replica_vs_original(G=base_graph, num_replicas=num_replicas, seed=1, params=params_MUSKETEER, metrics=metrics, title_infix='memoriless')
    

def paper_findings_one_metric(fname_data=None, alternative_algs=None):
    base_graph = graphutils.load_graph(path='data-epi/potterat.gml')
    base_graph.name = 'potterat'
    assert base_graph.number_of_nodes() == 250
    assert base_graph.number_of_edges() == 266
    num_replicas = 150
    #base_graph = graphutils.load_graph(path='data-social/911_unwted.gml')
    #base_graph.name = '911_krebs'
    
    print('Musketeer:')
    params  = {'verbose':False, 'edge_edit_rate':[0, 0, 0,   0, 0, 0.08], 'node_edit_rate':[0, 0, 0,   0, 0, 0.08], 'node_growth_rate':[0], 'dont_cutoff_leafs':True, 'enforce_connected':True, 'accept_chance_edges':0.8, 'new_edge_horizon':20, 'num_deletion_trials':20, 'skip_param_sanity_check':True}

    params['preserve_degree'] = True 
    params['epsilon']         = sum(params['edge_edit_rate'])+sum(params['node_edit_rate'])
    print('Params:')
    print(params)
    if base_graph.name == 'potterat':
        params['matrix'] = '0.629063, 0.584033; 0.668222, 0.13873'

    if alternative_algs==None:
        alternative_algs = []
        alternative_algs += [('MUSKETEER', algorithms.generate_graph)]
        alternative_algs += [('Edge swapping', alternatives.random_noise_replicate)]
        alternative_algs += [('Expected degree model', alternatives.expected_degree_replicate)]
        alternative_algs += [('Scale-free model', alternatives.scalefree_replicate)]
        alternative_algs += [('Kronecker model', alternatives.kronecker_replicate)]
    else:
        alternative_algs = [('MUSKETEER', algorithms.generate_graph)] + alternative_algs

    if fname_data != None:
        with open(fname_data, 'rb') as f:
            pkled_data = pickle.load(f)  
            replicas = pkled_data['replicas']
            params   = pkled_data['params']
    else:
        replicas = {}
        for alg_name,alg_func in alternative_algs:
            replicas[alg_name] = testers.replicate_graph(G=base_graph, generator_func=algorithms.generate_graph, num_replicas=num_replicas, params=params)

    #alternatives_metric(G=base_graph, replicas=replicas, 
    #                    metric={'name':'modularity', 'function':lambda G: community.louvain_modularity(G)})
    #alternatives_metric(G=base_graph, replicas=replicas, 
    #                    metric={'name':'total Deg*Deg', 'function':lambda G: nx.s_metric(G, normalized=False)})
    alternatives_metric(G=base_graph, replicas=replicas, 
                        metric={'name':'avg between central', 'function':lambda G: np.average(list(nx.betweenness_centrality(G, normalized=True).values()))})



def paper_findings_snapshots():
    base_graph = graphutils.load_graph(path='data-epi/potterat.gml')
    base_graph.name = 'potterat'
    assert base_graph.number_of_nodes() == 250
    assert base_graph.number_of_edges() == 266
    num_replicas = 150
    #base_graph = graphutils.load_graph(path='data-social/911_unwted.gml')
    #base_graph.name = '911_krebs'
    
    print('Musketeer:')
    params_MUSKETEER  = {'verbose':False, 'edge_edit_rate':[0.08, 0.07], 'node_edit_rate':[0.08, 0.07], 'node_growth_rate':[0], 'dont_cutoff_leafs':True, 'enforce_connected':True, 'accept_chance_edges':0.5, 'new_edge_horizon':20, 'num_deletion_trials':20, 'locality_bias_correction':[0, 0, 0, 0, 0]}


    ######################
    #     6th figure
    ######################
    metrics = graphutils.default_metrics[:]
    metrics  = [met for met in metrics if met['name'] not in ['avg flow closeness', 'avg eigvec centrality']]
    #metrics  = filter(lambda met: met['name'] not in ['avg harmonic shortest path'], metrics)

    replica_vs_original_snapshots(G=base_graph, num_snapshots=10, seed=1, params=params_MUSKETEER)
    


def paper_illustration():
    metrics = graphutils.default_metrics[:]
    metrics = [met for met in metrics if met['name'] not in ['avg flow closeness', 'powerlaw exp', 'density', 'num comps', 'harmonic mean path', 'mean ecc', 'average degree']]
    metrics = [met for met in metrics if met['name'] not in ['degree assortativity']]  #varies too much
    metrics = [met for met in metrics if met['optional'] == 0]

    ER = nx.erdos_renyi_graph(n=100, p=0.03)
    ER.name = 'ER'
    self_dissimilarity_analysis(G=ER, metrics=metrics)

    gpath     = 'output/%s'%ER.name + timeNow() + '.dot'
    gpath_fig = gpath[:-4] + '.eps'
    graphutils.write_graph(G=ER, path=gpath)
    visualizer_cmdl = 'sfdp -Nlabel="" -Nfontsize=0 -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Nfillcolor=/ylgnbu9/9 \
                            -Nstyle=filled -Nperipheries=0 -Gratio=0.75 -Gdpi=300 -Teps %s > %s &'%(gpath,gpath_fig)
    #http://www.graphviz.org/doc/info/colors.html
    print('Writing graph image: %s ..'%gpath_fig)
    retCode = os.system(visualizer_cmdl)

    #base_graph = graphutils.load_graph(path='data-epi/potterat.gml')
    #base_graph.name = 'potterat'
    #self_dissimilarity_analysis(G=base_graph, metrics=metrics)

    #WS = nx.watts_strogatz_graph(n=250, k=4, p=0.02)
    #WS.name = 'WS'
    #self_dissimilarity_analysis(G=WS, metrics=metrics)

    #BA = nx.barabasi_albert_graph(n=250, m=10, seed=42)
    #BA.name = 'BA'
    #self_dissimilarity_analysis(G=BA, metrics=metrics)

def paper_illustration_appendix():
    graphs = []
    graphs += [('sm-mesh', graphutils.load_graph(path='data/mesh33.gml'))]
    graphs += [('sm-watts', graphutils.load_graph(path='data-engineering/watts_strogatz98_power.elist'))]

    experiments = []
    experiments += [('Original network', 'orig',                            {'edge_edit_rate':[0],      'node_edit_rate':[0.],              'node_growth_rate':[0],})]
    experiments += [('Fine level changes', 'local',                         {'edge_edit_rate':[0.03],   'node_edit_rate':[0.],              'node_growth_rate':[0],})]
    experiments += [('Small number of coarse level changes', 'global_only', {'edge_edit_rate':[0, 0.02], 'node_edit_rate':[0],              'node_growth_rate':[0],})]
    experiments += [('Fine level changes with node growth', 'fine_growth',  {'edge_edit_rate':[0.03],   'node_edit_rate':[],                'node_growth_rate':[1.0],})]
    experiments += [('Several coarse level changes', 'several_global',      {'edge_edit_rate':[0.0, 0.06], 'node_edit_rate':[],             'node_growth_rate':[0],})]
    experiments += [('Coarse changes and node growth', 'coarse_growth',     {'edge_edit_rate':[0.0, 0.02], 'node_edit_rate':[0],            'node_growth_rate':[1.0],})]

    matplotlib.rc('text', usetex=False)
    
    print('Musketeer:')
    params  = {'verbose':False, 'dont_cutoff_leafs':True, 'enforce_connected':True, 'accept_chance_edges':1.0, 'new_edge_horizon':20, 'num_deletion_trials':20, 'locality_bias_correction':[0, 0, 0, 0, 0]}

    sfdp_default_cmd = 'sfdp -Goverlap="prism100" -Goverlap_scaling=-100 -Nlabel="" -Nwidth=0.01 -Nfixedsize=true -Nheight=0.01'
    
    for g_name, base_graph in graphs:
        print()
        gpath     = 'output/'+g_name+'_'#+timeNow()
        print(gpath)
        for exp_name, fname_suffix, param_update in experiments:
            print(exp_name)
            output_path = gpath+fname_suffix+'.dot'
            image_path  = output_path[:-4] + '.eps'
            stderr_path = output_path[:-4] + '.err.txt'
            
            params2 = params.copy()
            params2.update(param_update)
            replica = algorithms.generate_graph(original=base_graph, params=params2)
            graphutils.write_graph(G=replica, path=output_path)

            tmp_path = output_path+'_tmp'
            visualizer_cmdl = sfdp_default_cmd +' -Gdimen=2 -Txdot %s > %s 2> %s &'%(output_path,tmp_path,stderr_path)

            print('Writing graph with coordinates: %s ..'%tmp_path)
            sys.stdout.flush()
            retCode = os.system(visualizer_cmdl)

            print('Updating colors...')
            new_G = graphutils.update_visualization_attributes(nx.read_dot(tmp_path), verbose=False)
            print('Saving graph with updated layout and color: %s'%output_path)
            sys.stdout.flush()
            graphutils.write_graph(new_G, output_path)

            visualizer_cmdl = sfdp_default_cmd +' -Tpdf %s > %s 2> %s &'%(output_path,image_path,stderr_path)
            print('Writing graph image: %s ..'%image_path)
            sys.stdout.flush()
            retCode = os.system(visualizer_cmdl)

            #visualizer_cmdl = 'sfdp -Nlabel="" -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Nfillcolor=/ylgnbu9/9 -Nstyle=filled -Nperipheries=0 -Teps %s > %s &'%(output_path,image_path)
            #visualizer_cmdl = 'sfdp -Nlabel=""              -Nfixedsize=true               -Nfillcolor=/ylgnbu9/9 -Nstyle=filled -Nperipheries=0 -Teps %s > %s &'%(output_path,image_path)
            #http://www.graphviz.org/doc/info/colors.html
            #print 'Writing graph image: %s ..'%gpath_rep_fig
            #retCode = os.system(visualizer_cmdl)


def pickle_data(fname, data):
    try:
        os.mkdir('output')
    except:
        print('output/ cannot be created or already exists')

    with open(fname, 'wb') as outputFile:
        pickle.dump(data, outputFile)
        print() 
        print('Pickle: ' + fname + ' written!')
        return

    print('Unable to pickle...')

#seed=307269
#output/potterat_2013_05_15__18_11_46__base.dot.pdf 
#output/potterat_2013_05_15__18_10_51__node_edit_f.dot.pdf
#output/potterat_2013_05_15__18_11_26__node_edit_f.dot.pdf
#output/potterat_2013_05_15__18_11_46__node_edit_f.dot.pdf
def poster_figure_potterat(seed=None):
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)

    input_graph = graphutils.load_graph('data-epi/potterat.gml')
    G = nx.connected_component_subgraphs(input_graph)[0]
    G.name = 'potterat'
    for node in G:
        if '+' in G.node[node]['label']:
            G.node[node].update({'color':'red', 'style':'filled', 'label':''})  #'shape':'square', 
        elif '-' in G.node[node]['label']:
            G.node[node].update({'color':'black', 'style':'filled', 'label':''}) #'shape':'circle'}
        else:
            G.node[node].update({'color':'grey', 'style':'filled', 'label':''})
        G.node[node].pop('label')
        G.node[node].pop('shape')
        G.node[node].pop('x')
        G.node[node].pop('y')

    gpath = 'output/potterat_'+timeNow()

    params_default = {'maintain_node_attributes':True, 'algorithm':algorithms.musketeer_iterated_cycle, 'num_v_cycles':10,}
    p_sets = [{'name':'base',},
              {'name':'node_edit','new_edge_horizon':20, 'edge_edit_rate':[0.08,0.07  ], 'node_edit_rate':[0.08,0.07    ], 'node_growth_rate':[0       ],},
              #{'name':'node_edit_c','new_edge_horizon':20, 'edge_edit_rate':[0.0,  ], 'node_edit_rate':[0, 0, 0.05,     ], 'node_growth_rate':[0       ],},
              #{'name':'edge_edit_f','new_edge_horizon':20, 'edge_edit_rate':[0.1   ], 'node_edit_rate':[0.0,      ], 'node_growth_rate':[0       ],},
              #{'name':'node_grow','new_edge_horizon':20, 'edge_edit_rate':[0, 0, ], 'node_edit_rate':[0,            ], 'node_growth_rate':[0.2,     ],}, 
              #{'name':'edge_grow','new_edge_horizon':20, 'edge_edit_rate':[0, 0  ], 'node_edit_rate':[0             ], 'node_growth_rate':[        ], 'edge_growth_rate':[0.2],},  
              #{'name':'all_grow','new_edge_horizon':20, 'edge_edit_rate':[0, 0  ], 'node_edit_rate':[0             ], 'node_growth_rate':[0, 0.1        ], 'edge_growth_rate':[0, 0.1],},  
              ]#{'name':'kronecker', 'algorithm':alternatives.kronecker_replicate}]

    #params_default['algorithm'] = algorithms.musketeer_on_subgraphs
    for idx, p_set in enumerate(p_sets):
        params = params_default.copy()
        params.update(p_set)
        if 'name' in params:
            name = params.pop('name')
        else:
            name = str(idx)
        replica = algorithms.generate_graph(original=G, params=params)
        #replica = graphutils.color_new_nodes_and_edges(replica, G)

        print(replica.number_of_nodes())
        gpath_rep       = gpath + "__" + name + ".dot"
        gpath_rep_fig   = gpath_rep + ".pdf"
        graphutils.write_graph(G=replica, path=gpath_rep)
        #visualizer_cmdl = 'sfdp  -Nlabel="" -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Teps %s > %s &'
        print('Writing graph image: %s ..'%gpath_rep_fig)
        #visualizer_cmdl = 'sfdp -Goverlap="prism100" -Goverlap_scaling=-100 -Nlabel="" -Tpdf %s > %s '%(gpath_rep,gpath_rep_fig)
        visualizer_cmdl = 'sfdp -Nlabel="" -Nwidth=0.05 -Nheight=0.05 -Nfixedsize=true -Tpdf %s > %s '%(gpath_rep,gpath_rep_fig)
        retCode = os.system(visualizer_cmdl)
        time.sleep(0.3)
        subprocess.call(['xdg-open', gpath_rep_fig])



def presentation_figures():

    metrics_default = graphutils.default_metrics[:]
    metrics_default  = [met for met in metrics_default if met['name'] not in [#'avg flow closeness', 
    'avg eigvec centrality',
     'degree connectivity', 'harmonic mean path', 'mean ecc', 'powerlaw exp', 
     #'degree assortativity', 
    'avg between. central.', 
     ]]
    #metrics  = filter(lambda met: met['name'] not in ['avg harmonic shortest path'], metrics)

    #WS50 = nx.watts_strogatz_graph(50, 5, 0.1)
    #WS200= nx.watts_strogatz_graph(200, 5, 0.1)
    #Gpower = graphutils.load_graph('data-engineering/watts_strogatz98_power.elist')
    G_as = graphutils.load_graph('data-engineering/as-caida20071105.welist')
    Ggnutella = graphutils.load_graph('data-engineering/p2p-Gnutella09-gcc.elist')
    testers.replica_vs_original(G=Gpower, 
                                metrics=metrics_default,
            num_replicas=15, seed=21,  
            #params={'edge_edit_rate':[0.05, 0.05, 0.05, 0.05, 0.05, 0.05], 'node_edit_rate':[0.05, 0.05, 0.05, 0.05, 0.05, 0.05,], 'verbose':False, 'locality_bias_correction':[], })
            #params={'edge_edit_rate':[0, 0.10,  0, 0], 'node_edit_rate':[0.0, 0.10, 0, 0, 0.10], #WORKING
            #params={'edge_edit_rate':[0, 0.03, 0.07, 0], 'node_edit_rate':[0.0, 0, 0.7, 0, 0.10], 
            params={'edge_edit_rate':[0, 0, 0, 0.20, ], 'node_edit_rate':[0.0, 0, 0, 0.20, ],  
            #'algorithm':algorithms.musketeer_iterated_cycle, 'num_v_cycles':1, 
            'new_edge_horizon':20, 
            'verbose':False })


    '''
    base_graph = graphutils.load_graph('data-epi/adol.edgelist')
    base_graph.name = 'adolescent'
    gpath     = 'output/'+base_graph.name+'_'+timeNow()+'.dot'
    gpath_fig = gpath[:-3]+'png'
    graphutils.write_graph(G=base_graph, path=gpath)
    visualizer_cmdl = 'sfdp -Nlabel="" -Nfontsize=0 -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Nfillcolor=/ylgnbu9/9 \
                            -Nstyle=filled -Nperipheries=0 -Gratio=0.75 -Gdpi=300 -Tpng %s > %s &'%(gpath,gpath_fig)
    #http://www.graphviz.org/doc/info/colors.html
    print 'Writing graph image: %s ..'%gpath_fig

    retCode = os.system(visualizer_cmdl)

    #testers.replica_vs_original(G=graphutils.load_graph('data-social/911_unwted.gml'), metrics=metrics,
    #        num_replicas=50, seed=None,
    #        params={'edge_edit_rate':[0.08, 0.07], 'node_edit_rate':[0.08, 0.07], 'verbose':False, 'enforce_connected':True, 'locality_bias_correction':[0.9, 0.9], 'edge_welfare_fraction':1.0, 'edit_method':'alternating', })

    metrics = graphutils.default_metrics[:]
    metrics = filter(lambda met: met['name'] not in ['avg flow closeness', 'powerlaw exp'], metrics)

    testers.replica_vs_original(G=graphutils.load_graph('data-social/karate.adjlist'), metrics=metrics,
            num_replicas=50, seed=20,  
            params={'edge_edit_rate':[0.08, 0.07], 'edge_welfare_fraction':1.0, 'node_edit_rate':[0.08, 0.07], 'verbose':False, 'locality_bias_correction':[1.0, 1.0], })

    #testers.replica_vs_original(G=graphutils.load_graph('paperbenchmark/roadNet-TX.elist'), metrics=metrics,
    #        num_replicas=50, seed=None,
    #        params={'edge_edit_rate':[0.08, 0.07], 'node_edit_rate':[0.08, 0.07], 'verbose':False, 'enforce_connected':True, 'locality_bias_correction':[], 'edge_welfare_fraction':0.0, 'edit_method':'sequential', 'accept_chance_edges':0.0, 'deferential_detachment_factor':0})

    testers.replica_vs_original(G=graphutils.load_graph('paperbenchmark/yeast_cc1.elist'), metrics=metrics,
            num_replicas=50, seed=None,
            params={'edge_edit_rate':[0.08, 0.07], 'node_edit_rate':[0.08, 0.07], 'verbose':False, 'enforce_connected':True, 'locality_bias_correction':[], 'edge_welfare_fraction':0.0, 'edit_method':'alternating', 'accept_chance_edges':1.0, 'deferential_detachment_factor':1.0})


    #testers.replica_vs_original(G=graphutils.load_graph('data-cyber-small/gr2.gml'), num_replicas=150, seed=None, params={'edge_edit_rate':[0.08, 0.07], 'node_edit_rate':[0.08, 0.07], 'verbose':False, 'enforce_connected':True, 'locality_bias_correction':[], 'edge_welfare_fraction':0})
    
    #testers.replica_vs_original(G=graphutils.load_graph('data-social/ftp_core.elist'), num_replicas=30, seed=None, params={'edge_edit_rate':[0.08, 0.07], 'node_edit_rate':[0.08, 0.07], 'verbose':False, 'enforce_connected':True, 'locality_bias_correction':[], 'edit_method':'alternating', 'edge_welfare_fraction':0})
    '''

def presentation_figures_annotation(seed=1):
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)

    G = nx.grid_2d_graph(33,33)
    for node in G:
        if (node[0] + node[1]) % 2 == 0: #bipatite
            G.node[node] = {'color':'green', 'shape':'square', 'style':'filled'}
        else:# node[0] % 2 == 1:
            G.node[node] = {'color':'blue', 'shape':'circle'}
    graphutils.write_graph(G, 'data/mesh33-annotated.gml')
    gpath = 'output/mesh33an'

    params_default = {'maintain_node_attributes':True, 'algorithm':algorithms.musketeer_iterated_cycle, 'num_v_cycles':1,} #better use 10 cycles
    p_sets = [{'_name':'base',},
              #{'_name':'messy','edge_edit_rate':[0.0,  ], 'node_edit_rate':[0, 0, 0, 0.1,  ], 'node_growth_rate':[0       ], 'enforce_connected':False}, #MESS for testing of aggregates
              #{'_name':'node_edit_f','edge_edit_rate':[0.0,  ], 'node_edit_rate':[0.1,     ], 'node_growth_rate':[0       ],},
              #{'_name':'node_edit_c','edge_edit_rate':[0.0,  ], 'node_edit_rate':[0, 0, 0.05,     ], 'node_growth_rate':[0       ],},
              #{'_name':'edge_edit','edge_edit_rate':[0.1   ], 'node_edit_rate':[0.0,      ], 'node_growth_rate':[0       ],},
              #{'_name':'node_grow','edge_edit_rate':[0, 0, ], 'node_edit_rate':[0,            ], 'node_growth_rate':[0.2,     ],}, 
              #{'_name':'edge_grow','edge_edit_rate':[0, 0  ], 'node_edit_rate':[0             ], 'node_growth_rate':[        ], 'edge_growth_rate':[0.2],},  
              {'_name':'exp_deg', 'algorithm':alternatives.expected_degree_replicate}]
              #{'_name':'kronecker', 'algorithm':alternatives.kronecker_replicate}]

    #params_default['algorithm'] = algorithms.musketeer_on_subgraphs
    for idx, p_set in enumerate(p_sets):
        params = params_default.copy()
        params.update(p_set)
        replica = algorithms.generate_graph(original=G, params=params)
        print(replica.number_of_nodes())
        gpath_rep       = gpath + "__" + params.get('_name', str(idx)) + ".dot"
        gpath_rep_fig   = gpath_rep + ".pdf"
        graphutils.write_graph(G=replica, path=gpath_rep)
        #visualizer_cmdl = 'sfdp  -Nlabel="" -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Teps %s > %s &'
        print('Writing graph image: %s ..'%gpath_rep_fig)
        #visualizer_cmdl = 'sfdp -Goverlap="prism100" -Goverlap_scaling=-100 -Nlabel="" -Tpdf %s > %s '%(gpath_rep,gpath_rep_fig)
        visualizer_cmdl = 'sfdp -Nlabel="" -Nwidth=0.10 -Nheight=0.10 -Nfixedsize=true -Tpdf %s > %s '%(gpath_rep,gpath_rep_fig)
        retCode = os.system(visualizer_cmdl)
        time.sleep(0.3)
        subprocess.call(['xdg-open', gpath_rep_fig])


def presentation_figures_social(seed=1):
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)

    G = graphutils.load_graph('data-social/newman06_netscience.gml')
    base_graph = nx.connected_component_subgraphs(G)[0]
    #for node in G:
    #    if node[0] % 2 == 0:
    #        G.node[node] = {'color':'green', 'shape':'square', 'style':'filled'}
    #    else:# node[0] % 2 == 1:
    #        G.node[node] = {'color':'blue', 'shape':'circle'}
    gpath = 'output/netsci'

    params_default = {'maintain_node_attributes':True, 'algorithm':algorithms.musketeer_iterated_cycle, 'num_v_cycles':10,}
    p_sets = [{'name':'base',},
              #{'name':'node_edit_f','edge_edit_rate':[0.0,  ], 'node_edit_rate':[0.1,     ], 'node_growth_rate':[0       ],},
              #{'name':'node_edit_c','edge_edit_rate':[0.0,  ], 'node_edit_rate':[0, 0, 0.05,     ], 'node_growth_rate':[0       ],},
              #{'name':'edge_edit','edge_edit_rate':[0.1   ], 'node_edit_rate':[0.0,      ], 'node_growth_rate':[0       ],},
              #{'name':'node_grow','edge_edit_rate':[0, 0, ], 'node_edit_rate':[0,            ], 'node_growth_rate':[0.2,     ],}, 
              #{'name':'edge_grow','edge_edit_rate':[0, 0  ], 'node_edit_rate':[0             ], 'node_growth_rate':[        ], 'edge_growth_rate':[0.2],},  
              {'name':'all_grow','edge_edit_rate':[0, 0  ], 'node_edit_rate':[0             ], 'node_growth_rate':[0, 0.3        ], 'edge_growth_rate':[0, 0.1],},  
              ]#{'name':'kronecker', 'algorithm':alternatives.kronecker_replicate}]

    #params_default['algorithm'] = algorithms.musketeer_on_subgraphs
    for idx, p_set in enumerate(p_sets):
        params = params_default.copy()
        params.update(p_set)
        replica = algorithms.generate_graph(original=base_graph, params=params)
        replica = graphutils.color_new_nodes_and_edges(replica, base_graph)

        print(replica.number_of_nodes())
        gpath_rep       = gpath + "__" + params.get('name', str(idx)) + ".dot"
        gpath_rep_fig   = gpath_rep + ".pdf"
        graphutils.write_graph(G=replica, path=gpath_rep)
        #visualizer_cmdl = 'sfdp  -Nlabel="" -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Teps %s > %s &'
        print('Writing graph image: %s ..'%gpath_rep_fig)
        #visualizer_cmdl = 'sfdp -Goverlap="prism100" -Goverlap_scaling=-100 -Nlabel="" -Tpdf %s > %s '%(gpath_rep,gpath_rep_fig)
        visualizer_cmdl = 'sfdp -Nlabel="" -Nwidth=0.03 -Nheight=0.03 -Nfixedsize=true -Tpdf %s > %s '%(gpath_rep,gpath_rep_fig)
        retCode = os.system(visualizer_cmdl)
        time.sleep(0.3)
        subprocess.call(['xdg-open', gpath_rep_fig])


def presentation_figures_epi():
    base_graph = graphutils.load_graph('data-epi/proj.edgelist')
    base_graph = nx.connected_component_subgraphs(base_graph)[0]
    base_graph.name = 'proj90_gcc'

    gpath     = 'output/'+base_graph.name+'_'+timeNow()+'.dot'
    gpath_fig = gpath+'.pdf'
    graphutils.write_graph(G=base_graph, path=gpath)
    print('Writing graph image: %s ..'%gpath_fig)
    visualizer_cmdl = 'sfdp -Nlabel="" -Nwidth=0.03 -Nheight=0.03 -Nfixedsize=true -Nfillcolor=/ylgnbu9/9 -Nshape=point -Tpdf %s > %s &'%(gpath,gpath_fig)
    retCode = os.system(visualizer_cmdl)
    #http://www.graphviz.org/doc/info/colors.html

    params_default  = {'verbose':False, 
            'dont_cutoff_leafs':True, #NOTE: this is important
            'new_edge_horizon':10, 'num_deletion_trials':20, 'locality_bias_correction':[0,],
            'accept_chance_edges':0.5, 
            'enforce_connected':True, 
            #'enforce_connected':False, 
            'component_is_edited': [True] + 1000*[False],
            #'deferential_detachment_factor':0.0, #1.0,
            }

    p_sets = [{'edge_edit_rate':[0.04,     ], 'node_edit_rate':[0.04,     ], 'node_growth_rate':[0],},
              {'edge_edit_rate':[0, 0, 0.04], 'node_edit_rate':[0, 0, 0.04], 'node_growth_rate':[0],}, 
              {'edge_edit_rate':[0, 0      ], 'node_edit_rate':[0         ], 'node_growth_rate':[0.1],},  ]

    #params_default['algorithm'] = algorithms.musketeer_on_subgraphs
    for idx, p_set in enumerate(p_sets):
        params = params_default.copy()
        params.update(p_set)
        replica = algorithms.generate_graph(original=base_graph, params=params)
        print(replica.number_of_nodes())
        gpath_rep = gpath + "__" + str(idx) + ".dot"
        gpath_rep_fig = gpath_rep + ".pdf"
        graphutils.write_graph(G=replica, path=gpath_rep)
        #visualizer_cmdl = 'sfdp  -Nlabel="" -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Teps %s > %s &'
        print('Writing graph image: %s ..'%gpath_rep_fig)
        visualizer_cmdl = 'sfdp -Nlabel="" -Nwidth=0.03 -Nheight=0.03 -Nfixedsize=true -Nfillcolor=/ylgnbu9/9 -Nshape=point -Tpdf %s > %s &'%(gpath_rep,gpath_rep_fig)
        #visualizer_cmdl = 'sfdp -Nlabel="" -Nfontsize=0 -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Ncolor=blue \
        #                   -Nshape=point  -Npenwidth=0 -Nperipheries=0 -Gratio=0.75 -Gdpi=300 -Tpdf %s > %s &'%(gpath_rep,gpath_rep_fig)
        retCode = os.system(visualizer_cmdl)




def replica_vs_original_degree_distribution(seed=None, figpath=None, show_replicas=False, generator_func=None, G=None, params=None, num_replicas = 150):
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)

    if generator_func==None:
        generator_func=algorithms.generate_graph

    if G==None:
        G = nx.barabasi_albert_graph(n=200, m=10)
        G.name = 'Barabasi-Albert n=200, m=10'
        #graphs.append(graphutils.load_graph('data-social/911_unwted.gml'))
        #graphs[-1].name = '911 Network'
        #graphs.append(graphutils.load_graph('data-cyber-small/gr2.gml'))
        #graphs[-1].name = 'Cyber GR2'
        #graphs.append(graphutils.load_graph('data-cyber-small/gr1.gml'))
        #graphs[-1].name = 'Cyber GR1'

    num_replicas =  150

    if params == None:
        params  = {'verbose':False, 'edge_edit_rate':[0.05, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01], 'node_edit_rate':[0.05, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01], 'node_growth_rate':[0., 0], 'locality_bias_correction':[0.], 'enforce_connected':True, 'accept_chance_edges':1.0}
    print('Params:')
    print(params)
    vals_of_replicas = []
    #nums_of_nodes    = []
    vals_of_graph    = []
    for replica_idx in range(num_replicas):
        replica     = generator_func(G, params=params)
        num_nodes   = float(replica.number_of_nodes())
        vals_of_replicas.append([v/num_nodes for v in nx.degree_histogram(replica)])
        sys.stdout.write('.')
        #nums_of_nodes.append(num_nodes)
    pylab.figure()
    pylab.hold(True)
    max_dist_support = max(len(vals) for vals in vals_of_replicas)
    for i, vals in enumerate(vals_of_replicas):
        #pylab.plot(vals, '-.', color='0.5', markersize=1)
        vals += (max_dist_support-len(vals))*[0.]
        vals_of_replicas[i] = 1 - np.cumsum(vals) 

    original_num_nodes = G.number_of_nodes()
    original_histogram = np.array(nx.degree_histogram(G))/float(original_num_nodes)
    xvals = 1 - np.cumsum(original_histogram)
    pylab.plot(xvals, 'o-', color='b', linewidth=2., label='empirical')

    vals_of_replicas = np.array(vals_of_replicas)
    mean_of_replicas = vals_of_replicas.mean(axis=0)
    std_of_replicas  = vals_of_replicas.std(axis=0, ddof=1)
    
    stepping = 1
    pylab.plot(list(range(max_dist_support)), mean_of_replicas, '-', color='g', linewidth=1., label='generated networks')
    pylab.errorbar(list(range(max_dist_support))[::stepping], mean_of_replicas[::stepping], yerr=std_of_replicas[::stepping], color='r', linestyle='', ecolor='g', capsize=4)[0]
    
    pylab.ylabel(r'$\bf{P}[$Degree $> k]$', fontsize='20')
    #pylab.ylabel(r'$\mathbb{P}[Degree > k]$', fontsize='20')
    pylab.xlabel(r'$k$', fontsize='20')#, x=0.1)
    #pylab.ylabel('Normed by mean value', rotation=90, fontsize='20')
    #pylab.title(G.name)
    pylab.legend(loc='best')
    pylab.ylim(-0.01,1.01)

    if figpath == None:
        figpath = 'output/replica_vs_original_degree_hist_'+G.name+'_'+str(seed)+'__'+timeNow()
        figpath = clean_path(figpath)
    save_figure_helper(figpath)
    pylab.hold(False)



def replica_vs_original_snapshots(G, params, verbose=True, num_snapshots=10, seed=None, metrics=None, figure=None, norms=None, markersize=7):
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)

    snapshot_params = params.copy()
    snapshot_params = {'edge_edit_rate':np.array(params['edge_edit_rate'])/num_snapshots,
                       'node_edit_rate':np.array(params['node_edit_rate'])/num_snapshots,
                       'node_growth_rate':np.array(params['node_growth_rate'])/num_snapshots,
                       }

    snapshot_data = {}
    clean_names = {'num nodes': 'num nodes', 'num edges':'num edges', 'clustering':'clustering', 'average degree':'avg degree', 'degree assortativity':'degree\nassortativity',  'degree connectivity':'degree\nconnectivity', 'total deg*deg':'total deg*deg\nassortativity', 's-metric':'s-metric', 'mean ecc':'avg\neccentricity', 'num comps':'num comps', 'L eigenvalue sum':'L eigen-\nvalue sum', 'average shortest path':'avg distance', 'harmonic mean path':'harmonic avg distance', 'avg flow closeness':'avg flow\ncloseness', 'avg eigvec centrality':'avg eigenvec. centrality', 'avg between. central.':'avg between. centrality', 'modularity':'modularity'}
    if metrics == None:
        metrics = [met for met in graphutils.default_metrics[:] if met['name'] not in ['avg flow closeness', 'powerlaw exp', 'density', 'num comps', 'mean ecc',]]
        metrics = [met for met in metrics if met['name'] not in ['degree assortativity']]  #varies too much
        metrics = [met for met in metrics if met['name'] not in ['degree connectivity']]  #varies too much
        metrics = [met for met in metrics if met['optional'] == 0]
    for met_info in metrics:
        snapshot_data[met_info['name']] = []

    params['num_snapshots'] = num_snapshots
    final_graph = algorithms.musketeer_snapshots(original=G, params=snapshot_params)
    graphs = final_graph.snapshots

    for snapshot in graphs:
        #print 'nn: %d, edges: %d'%(G_coarse.number_of_nodes(),G_coarse.number_of_edges())
        sys.stdout.write('.')
        for met_info in metrics:
            met_name = met_info['name']
            met_func = met_info['function']
            try:
                if met_name in ['harmonic mean path', 'average shortest path']:
                    snapshot_data[met_name] += [float(met_func(snapshot, max_num_sources=300))]
                else:
                    snapshot_data[met_name] += [float(met_func(snapshot))]
            except:
                continue
    num_levels = len(list(snapshot_data.values())[0])

    if norms ==None:
        norms = {}
        for idx,metric in enumerate(metrics):
            met_name = metric['name']
            try:
                if met_name in ['harmonic mean path', 'average shortest path']:
                    norms[met_name] = metric['function'](G, max_num_sources=300)
                else:
                    norms[met_name] = metric['function'](G)
                #norms[met_name] = np.average(snapshot_data[met_name])
            except Exception as inst:
                print('Warning: could not plot '+met_name + ': '+str(inst))
    if figure==None:
        pylab.figure()
        pylab.hold(True)    
    #curvetypes = ['r-', 'g-', 'b-', 'k-', 'c-.', 'y-', 'm-D', 'b-.', 'r-.', 'c-']
    curvetypes = ['ro-', 'g^-', 'bv-', 'ks-', 'cp-', 'yH-', 'mD-', 'b3-', 'r4-', 'c1-']
    linewidths = np.arange(1, 2, 0.1)
    for idx,metric in enumerate(metrics):
        met_name = metric['name']
        try:
            data = snapshot_data[met_name]
            pylab.plot(np.arange(len(data)),# + npr.rand(len(data))/20.,  #jitter
                       np.array(data)/norms[met_name], 
                       curvetypes[idx%len(curvetypes)], linewidth=linewidths[idx%len(linewidths)], label=clean_names.get(met_name, met_name), markersize=markersize)
        except Exception as inst:
            print('Warning: could not plot '+met_name + ': '+str(inst))
            raise 
    if figure==None:
        pylab.ylim(0.45, 1.55)
        pylab.vlines(list(range(1,10)), ymin=-0.05, ymax=2., linestyles='dashed')
        pylab.xlim(-0.05, num_snapshots+0.05)
        pylab.xlabel('Time (steps)', fontsize='20')
        #pylab.ylabel('Normalized value', fontsize='20')
        #pylab.title(G.name)
        pylab.legend(prop={'size':12}, ncol=3, loc='best')
        pylab.xticks(list(range(num_levels)),list(range(num_levels)))
        fname = 'output/replica_dynamics_'+getattr(G, 'name', '')+timeNow()
        fname = clean_path(fname)
        save_figure_helper(fname)

    #pylab.show()

    return snapshot_data, norms


def running_time_analysis(verbose=True):
    print()
    print('running_time_analysis')
    print()

    base_density = 0.01
    num_trials = 5
    alg_params = {}
    alg_params['verbose'] = False
    alg_params['edge_edit_rate'] = [0.05 for i in range(100)]
    alg_params['node_edit_rate'] = [0.05 for i in range(100)]
    alg_params['node_growth_rate'] = [0.05 for i in range(100)]
 
    num_nodes_array = [int(10**i) for i in np.arange(2, 4, 0.5)]
    density_array   = [base_density*i for i in range(1, 5)]

    solns   = [] 
    for i,nn in enumerate(num_nodes_array):
        solns.append([])
        for j,dens in enumerate(density_array):
            trial_outcomes = []
            for trial in range(num_trials):
                G = nx.erdos_renyi_graph(n=nn, p=dens)
                t0 = time.time()
                replica = algorithms.generate_graph(original=G, params=alg_params)
                solve_time = time.time() - t0
                trial_outcomes.append(solve_time)
            solns[i].append(np.average(trial_outcomes))
            print('nn=%d,p=%.3f:\t avg t=%.2f'%(nn,dens,solns[i][-1]))

    pickle_data(fname='output/runtime__vs_num_nodes' + timeNow() + '.pkl', data={'solns':solns, 'num_nodes_array':num_nodes_array, 'density_array':density_array})

    pylab.figure()
    pylab.hold(True)    
    curvetypes = ['r-', 'g-', 'b-', 'k-', 'c-.', 'y-', 'm-D', 'b-.', 'r-.']
    linewidths = [4.0, 2.5, 1.0, 0.5, 2.0, 2.0, 1.0, 5.0]
    for idx,nn in enumerate(num_nodes_array):
        pylab.plot(density_array, solns[idx], curvetypes[idx], linewidth=linewidths[idx], label='nodes=%d'%nn)
    pylab.xlabel('density')
    pylab.ylabel('running time (seconds)')
    pylab.title('')
    pylab.legend(loc='best')
    #pylab.xticks(range(num_levels),range(num_levels))
    pylab.savefig('output/running_time_ER_'+timeNow())

    #with open('output/runtime__vs_num_nodes2012_01_03__07_44_03.pkl', 'rb') as f:  ret = cPickle.load(f)
    #compute the slope by solving least square on log-log
    #see example in np.linalg.lstsq(a, b, rcond=-1)
    #suppose y = c0 + c1 * (x**m)
    c00 = solns[0][0]
    c01 = solns[0][-1]
    x = np.log(num_nodes_array[1:])
    y1 = np.log([solns[i][0]  - c00 for i in range(1, len(solns))])  #lowest density
    y2 = np.log([solns[i][-1] - c01 for i in range(1, len(solns))]) #highest density
    A = np.vstack([x, np.ones(len(x))]).T
    slope_term1, const_term1 = np.linalg.lstsq(A, y1)[0]
    slope_term2, const_term2 = np.linalg.lstsq(A, y2)[0]
    print('best fit exponent: %.4f to %.4f'%(slope_term1,slope_term2))
    #pylab.show()

    return solns


def save_figure_helper(fpath):
    if matplotlib.get_backend() == 'pdf':
        pylab.savefig(fpath + '.pdf')
        print('Written: %s'%fpath + '.pdf')
        print('Converting to eps...')
        os.system('pdftops -eps ' + fpath+'.pdf')
        os.rename(fpath+'.eps', fpath+'_.eps') #missing bounding box
        os.system('eps2eps ' + fpath+'_.eps' + '  ' + fpath+'.eps')
        os.remove(fpath+'_.eps')
    elif matplotlib.get_backend() == 'ps':
        pylab.savefig(fpath + '.eps')
        print('Written: %s'%fpath + '.eps')
        print('Converting to pdf...')
        os.system('epstopdf ' + fpath+'.eps')
    else:
        print('Saving of figure to unknown extension.  Backend: %s'%matplotlib.get_backend())
        pylab.savefig(fpath)
        print('Written: %s'%fpath + 'SOMETHING')


def self_dissimilarity_analysis(G, verbose=True, metrics=None, figure=None, norms=None, markersize=7):
    '''
    show how various graph metrics change with coarsening
    '''

    G_coarse = G.copy()
    params   = {}
    scaling_data = {}
    if metrics == None:
        for met_info in graphutils.default_metrics:
            if met_info['optional'] > 0:
                continue
    for met_info in metrics:
        scaling_data[met_info['name']] = []

    while G_coarse.number_of_nodes() > 1 and G_coarse.number_of_edges() > 1 :
        #print 'nn: %d, edges: %d'%(G_coarse.number_of_nodes(),G_coarse.number_of_edges())
        for met_info in metrics:
            met_name = met_info['name']
            met_func = met_info['function']
            try:
                scaling_data[met_name] += [float(met_func(G_coarse))]
            except:
                continue
        G_coarse = algorithms.do_coarsen(G_coarse, params)[0]
    num_levels = len(list(scaling_data.values())[0])

    if norms ==None:
        norms = {}
        for idx,metric in enumerate(metrics):
            met_name = metric['name']
            try:
                norms[met_name] = max(scaling_data[met_name])
            except Exception as inst:
                print('Warning: could not plot '+met_name + ': '+str(inst))
    if figure==None:
        pylab.figure()
        pylab.hold(True)    
    #curvetypes = ['r-', 'g-', 'b-', 'k-', 'c-.', 'y-', 'm-D', 'b-.', 'r-.', 'c-']
    curvetypes = ['ro--', 'g^--', 'bv-.', 'ks--', 'cp--', 'yH-.', 'mD-.', 'b3-.', 'r4-.', 'c1-']
    linewidths = np.arange(1, 2, 0.1)
    for idx,metric in enumerate(metrics):
        met_name = metric['name']
        try:
            data = scaling_data[met_name]
            pylab.plot(np.arange(len(data)) + npr.rand(len(data))/20.,  #jitter
                       np.array(data)/norms[met_name], 
                       curvetypes[idx%len(curvetypes)], linewidth=linewidths[idx%len(linewidths)], label=met_name, markersize=markersize)
        except Exception as inst:
            print('Warning: could not plot '+met_name + ': '+str(inst))
            raise 
    if figure==None:
        pylab.xlabel('scale (finest to coarsest)', fontsize='20')
        #pylab.ylabel('Normalized value', fontsize='20')
        #pylab.title(G.name)
        pylab.legend(loc='best')
        pylab.xticks(list(range(num_levels)),list(range(num_levels)))
        pylab.savefig('output/self_asymmetry_'+getattr(G, 'name', '')+timeNow()+'.eps')

    #pylab.show()

    return scaling_data, norms


def self_dissimilarity_demo():    
    WS200 = nx.watts_strogatz_graph(n=200, k=4, p=0.02)
    self_dissimilarity_analysis(G=WS200)

    ER200 = nx.erdos_renyi_graph(n=200, p=0.05)
    self_dissimilarity_analysis(G=ER200)

    BA200 = nx.barabasi_albert_graph(n=200, m=10, seed=42)
    self_dissimilarity_analysis(G=BA200)


def sims_on_powergrids(seed=1):
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)

    num_replicas = 10
    
    params_default  = {'verbose':False, 
                        'edge_edit_rate':[0.0], 'node_edit_rate':[0, 0, 0, 0.10],
                        'enforce_connected':True,
                        'num_v_cycles':5,
                      }
    params_default['algorithm'] = (algorithms.musketeer_on_subgraphs, algorithms.musketeer_iterated_cycle,)
    #params_for_big = params_default.copy()
    #params_for_big['num_v_cycles'] = 5
    #params_for_big['algorithm'] = algorithms.musketeer_on_subgraphs
    #metrics_for_big = [k for k in graphutils.default_metrics if k['name'] in ['num nodes', 'num edges', 'clustering', 'total deg*deg', 'harmonic mean path']]

    metrics_default = graphutils.default_metrics[:]
    metrics_default  = [met for met in metrics_default if met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree connectivity', 'degree assortativity',  'average shortest path', 'mean ecc', 'powerlaw exp', ]]
    #metrics  = filter(lambda met: met['name'] not in ['avg harmonic shortest path'], metrics)

    problems = [
                {'graph_path':'data-engineering/polish_power/Polish110Volts.edges', 'graph_params':params_default, 'graph_metrics':metrics_default,},
                #{'graph_path':'data-engineering/polish_power/Polish220Volts.edges', 'graph_params':params_default, 'graph_metrics':metrics_default,},
                #{'graph_path':'data-engineering/polish_power/Polish400Volts.edges', 'graph_params':params_default, 'graph_metrics':metrics_default,},
                ]

    for problem_info in problems:
        graph_path      = problem_info['graph_path']
        graph_params    = problem_info.get('graph_params', params_default)
        graph_metrics   = problem_info.get('graph_metrics', metrics_default)
        graph_num_replicas = problem_info.get('num_replicas', num_replicas)
        base_graph      = graphutils.load_graph(path=graph_path, params={'verbose':False})
        print('Size nn: %d, ne: %d'%(nx.number_of_nodes(base_graph), nx.number_of_edges(base_graph)))
        base_graph = nx.Graph(base_graph)
        print(graph_path, nx.number_of_nodes(base_graph), nx.number_of_edges(base_graph))
        
        base_graph.name = graph_path.split('/')[-1]

        gpath     = 'output/'+base_graph.name+'_'+timeNow()+'.dot'
        gpath_fig = gpath[:-3]+'eps'
        graphutils.write_graph(G=base_graph, path=gpath)
        try:
            visualizer_cmdl = 'sfdp  -Nlabel="" -Nwidth=0.03 -Nfixedsize=true -Nheight=0.03 -Teps %s > %s &'%(gpath,gpath_fig)
            print('Writing graph image: %s ..'%gpath_fig)
            retCode = os.system(visualizer_cmdl)
        except Exception as inst:
            print('Warning: could not plot '+graph_path + ': '+str(inst))   
            exc_traceback = sys.exc_info()[2]
            print(str(inst) + "\n" + str(traceback.format_tb(exc_traceback)).replace('\\n', '\n'))
    
        try:
            replicas = testers.replicate_graph(G=base_graph, generator_func=algorithms.generate_graph, num_replicas=num_replicas, params=graph_params)
            for rep_num, rep in enumerate(replicas):
                rep_path = 'output/'+base_graph.name+'_rep=%d'%(rep_num)+timeNow()+'.edges'
                rep = nx.convert_node_labels_to_integers(rep) 
                #graphutils.write_graph(G=rep, path=rep_path)
                nx.write_edgelist(G=rep, path=rep_path, data=False)
                print('Written: ' + rep_path)

            #testers.replica_vs_original(G=base_graph, num_replicas=graph_num_replicas, seed=seed, params=param_set, metrics=graph_metrics, title_infix='', n_jobs=-1)
        except Exception as inst:
            print('Warning: could not generate or evaluate '+graph_path + ': '+str(inst))
            exc_traceback = sys.exc_info()[2]
            print(str(inst) + "\n" + str(traceback.format_tb(exc_traceback)).replace('\\n', '\n'))


def test1():
    G = nx.connected_component_subgraphs(nx.erdos_renyi_graph(100, 0.2))[0]
    #G = nx.watts_strogatz_graph(20,4, 0.1)
    #params = {'error_rate':G.number_of_nodes()*[0.]}
    #params = {'error_rate':[0.03,] + 100*[0.]}
    params ={}
    new_G = generate_graph(G, params=params)
    
    pylab.figure()
    nx.draw(G)
    pylab.figure()
    nx.draw(new_G)
    benchmarks.compare_nets(G, new_G, params={'verbose':True})
    #pylab.show()
    
def test2():
    #on a 2D grid, artificially define seeds and run uncoarsen several times: it should restore the grid
    pass


def test3():
    pass #under a singular error rate (all focused on a particular level) moving the error up the levels should increase change in the graph

def testGrid():
    #Musketeer is not the right algorithm for grids, but we can see the effect of noise more clearly here
    G = nx.grid_2d_graph(10,5)
    pos = {}
    for node in G:
        pos[node] = (node[0] + npr.rand()/15., node[1] + npr.rand()/15.)  #jitter the positions
    
    #params = {'error_rate':[0.25,] + 100*[0.], 'verbose':True, 'debug_info':True}
    #params = {'error_rate':[0.,0., 0.25] + 100*[0.]}
    #params = {'error_rate':[0.,0.,0.05] + 100*[0.]}
    params = {'error_rate':100*[0.2]}
    #params ={}

    new_G = generate_graph(G, params=params)

    pylab.figure()
    nx.draw(new_G, pos=pos)
    benchmarks.compare_nets(G, new_G, params={'verbose':True})
    #pylab.show()

def find_differences(G,G_new):
    deleted_edges = 0
    new_edges = 0
    for edge in G.edges():
        if G_new.has_edge(edge[0],edge[1]) or G_new.has_edge(edge[1],edge[0]):
            continue
        else:
            deleted_edges+=1
    for edge in G_new.edges():
        if G.has_edge(edge[0],edge[1]) or G.has_edge(edge[1],edge[0]):
            continue
        else:
            new_edges+=1
    return new_edges,deleted_edges;

def test_planar_1():
    num_replicas = 30
    #params_MUSKETEER = {'verbose': False, '-k':True,'edge_edit_rate': [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01],
     #                   'node_growth_rate': [0],
      #                  'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
       #                 'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_MUSKETEER = {'verbose': False, '-k': True,
                        'edge_edit_rate': [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05],
                        'node_growth_rate': [0],
                        'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                        'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_coarse = {'verbose': False, '-k': True,
                        'edge_edit_rate': [ 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.10],
                        'node_growth_rate': [0],
                        'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                        'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_fine = {'verbose': False, '-k': True,
                     'edge_edit_rate': [0.1, 0.1, 0.10,0.0, 0.0, 0.0,],
                     'node_growth_rate': [0],
                     'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                     'num_trial_particles': 50, 'num_insertion_trials': 50}
    metrics = graphutils.default_metrics[:]
    metrics = [met for met in metrics if
               met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree assortativity',
                                   'degree connectivity']]

    ######################
    #     3rd figure
    ######################
    print('Barabasi-Albert')
    params_BA = params_MUSKETEER.copy()
    params_BA['locality_bias_correction'] = []
    #metrics = [met for met in graphutils.default_metrics if met['name'] == 'num edges']
    metrics[-1]['optional'] = 0
    #BA_graph = nx.barabasi_albert_graph(n=20, m=3, seed=42)
    BA_graph = nx.read_edgelist('/home/varsha/Documents/final_results/input/sample_graph_1.edgelist')
    #n=20
    #w = {i: random.expovariate(5.0) for i in range(n)}
    #BA_graph = nx.geographical_threshold_graph(20, 50, weight=w)

    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_BA, metrics=metrics,graph_name = 'sample_graph_1', generator = "musketeer_all",n_jobs=10)
    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_coarse, metrics=metrics,
                                graph_name='sample_graph_1', generator="musketeer_coarse", n_jobs=10)
    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_fine, metrics=metrics,
                                graph_name='sample_graph_1', generator="musketeer_fine", n_jobs=10)

def test_planar_road_network():
    num_replicas = 30
    #params_MUSKETEER = {'verbose': False, '-k':True,'edge_edit_rate': [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01],
     #                   'node_growth_rate': [0],
      #                  'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
       #                 'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_MUSKETEER = {'verbose': False, '-k': True,
                        'edge_edit_rate': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0],
                        'node_growth_rate': [0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15],
                        'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                        'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_coarse = {'verbose': False, '-k': True,
                        'edge_edit_rate': [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0],
                        'node_growth_rate': [0.0, 0.0, 0.0,0.0, 0.30, 0.30, 0.30, 0.30],
                        'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                        'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_fine = {'verbose': False, '-k': True,
                     'edge_edit_rate': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0],
                     'node_growth_rate': [0.30, 0.30, 0.30, 0.30,0.0, 0.0, 0.0,0.0],
                     'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                     'num_trial_particles': 50, 'num_insertion_trials': 50}
    metrics = graphutils.default_metrics[:]
    metrics = [met for met in metrics if
               met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree assortativity',
                                   'degree connectivity']]

    ######################
    #     3rd figure
    ######################
    print('Barabasi-Albert')
    params_BA = params_MUSKETEER.copy()
    params_BA['locality_bias_correction'] = []
    #metrics = [met for met in graphutils.default_metrics if met['name'] == 'num edges']
    metrics[-1]['optional'] = 0
    #BA_graph = nx.barabasi_albert_graph(n=20, m=3, seed=42)
    BA_graph = nx.read_edgelist('/home/varsha/Documents/final_results/input/sample_graph_1.edgelist')
    #n=20
    #w = {i: random.expovariate(5.0) for i in range(n)}
    #BA_graph = nx.geographical_threshold_graph(20, 50, weight=w)

    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_BA,
                                metrics=metrics,graph_name = 'sample_graph_1', generator = "musketeer_all",n_jobs=10)
    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_coarse, metrics=metrics,
                                graph_name='sample_graph_1', generator="musketeer_coarse", n_jobs=10)
    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_fine, metrics=metrics,
                                graph_name='sample_graph_1', generator="musketeer_fine", n_jobs=10)

def test_planar_nasa_planar():
    num_replicas = 30
    #params_MUSKETEER = {'verbose': False, '-k':True,'edge_edit_rate': [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01],
     #                   'node_growth_rate': [0],
      #                  'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
       #                 'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_MUSKETEER = {'verbose': False, '-k': True,
                        'edge_edit_rate': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0],
                        'node_growth_rate': [0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15,0.15, 0.15],
                        'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                        'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_coarse = {'verbose': False, '-k': True,
                        'edge_edit_rate': [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0],
                        'node_growth_rate': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.30, 0.30, 0.30, 0.30],
                        'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                        'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_fine = {'verbose': False, '-k': True,
                     'edge_edit_rate': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0],
                     'node_growth_rate': [0.30, 0.30, 0.30, 0.30,0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                     'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                     'num_trial_particles': 50, 'num_insertion_trials': 50}
    metrics = graphutils.default_metrics[:]
    metrics = [met for met in metrics if
               met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree assortativity',
                                   'degree connectivity']]

    ######################
    #     3rd figure
    ######################
    print('Barabasi-Albert')
    params_BA = params_MUSKETEER.copy()
    params_BA['locality_bias_correction'] = []
    #metrics = [met for met in graphutils.default_metrics if met['name'] == 'num edges']
    metrics[-1]['optional'] = 0
    #BA_graph = nx.barabasi_albert_graph(n=20, m=3, seed=42)
    BA_graph = nx.read_edgelist('/home/varsha/Documents/final_results/input/nasa_planar.edges')
    #n=20
    #w = {i: random.expovariate(5.0) for i in range(n)}
    #BA_graph = nx.geographical_threshold_graph(20, 50, weight=w)

    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_BA, metrics=metrics,
                                graph_name = 'nasa_planar', generator = "musketeer_all_rescale",n_jobs=5)
    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_coarse, metrics=metrics,
                                graph_name='nasa_planar', generator="musketeer_coarse_rescale", n_jobs=5)
    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_fine, metrics=metrics,
                                graph_name='nasa_planar', generator="musketeer_fine_rescale", n_jobs=5)

def test_planar_powergrid():
    num_replicas = 30
    #params_MUSKETEER = {'verbose': False, '-k':True,'edge_edit_rate': [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01],
     #                   'node_growth_rate': [0],
      #                  'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
       #                 'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_MUSKETEER = {'verbose': False, '-k': True,
                        'edge_edit_rate': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0],
                        'node_growth_rate': [0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15,0.15],
                        'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                        'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_coarse = {'verbose': False, '-k': True,
                        'edge_edit_rate': [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0],
                        'node_growth_rate': [ 0.0, 0.0, 0.0, 0.0, 0.30, 0.30, 0.30, 0.30],
                        'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                        'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_fine = {'verbose': False, '-k': True,
                     'edge_edit_rate': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0],
                     'node_growth_rate': [0.30, 0.30, 0.30, 0.30, 0.0, 0.0, 0.0, 0.0],
                     'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                     'num_trial_particles': 50, 'num_insertion_trials': 50}
    metrics = graphutils.default_metrics[:]
    metrics = [met for met in metrics if
               met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree assortativity',
                                   'degree connectivity']]

    ######################
    #     3rd figure
    ######################
    print('Barabasi-Albert')
    params_BA = params_MUSKETEER.copy()
    params_BA['locality_bias_correction'] = []
    #metrics = [met for met in graphutils.default_metrics if met['name'] == 'num edges']
    metrics[-1]['optional'] = 0
    #BA_graph = nx.barabasi_albert_graph(n=20, m=3, seed=42)
    BA_graph = nx.read_edgelist('/home/varsha/Documents/final_results/input/powergrid_subgraph.edges')
    #n=20
    #w = {i: random.expovariate(5.0) for i in range(n)}
    #BA_graph = nx.geographical_threshold_graph(20, 50, weight=w)

    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_BA, metrics=metrics,
                                graph_name = 'powergrid', generator = "musketeer_all_rescale",n_jobs=5)
    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_coarse, metrics=metrics,
                                graph_name='powergrid', generator="musketeer_coarse_rescale", n_jobs=5)
    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_fine, metrics=metrics,
                                graph_name='powergrid', generator="musketeer_fine_rescale", n_jobs=5)

def test_planar_water():
    num_replicas = 30
    #params_MUSKETEER = {'verbose': False, '-k':True,'edge_edit_rate': [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01],
     #                   'node_growth_rate': [0],
      #                  'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
       #                 'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_MUSKETEER = {'verbose': False, '-k': True,
                        'edge_edit_rate': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        'node_growth_rate': [0.15, 0.15, 0.15, 0.15, 0.15, 0.15,0.15],
                        'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                        'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_coarse = {'verbose': False, '-k': True,
                        'edge_edit_rate': [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        'node_growth_rate': [0.0, 0.0,0.0, 0.3, 0.3, 0.30,.30],
                        'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                        'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_fine = {'verbose': False, '-k': True,
                     'edge_edit_rate': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                     'node_growth_rate': [0.30, 0.30,0.30,0.30, 0.0, 0.0, 0.0],
                     'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                     'num_trial_particles': 50, 'num_insertion_trials': 50}
    metrics = graphutils.default_metrics[:]
    metrics = [met for met in metrics if
               met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree assortativity',
                                   'degree connectivity']]

    ######################
    #     3rd figure
    ######################
    print('Barabasi-Albert')
    params_BA = params_MUSKETEER.copy()
    params_BA['locality_bias_correction'] = []
    #metrics = [met for met in graphutils.default_metrics if met['name'] == 'num edges']
    metrics[-1]['optional'] = 0
    #BA_graph = nx.barabasi_albert_graph(n=20, m=3, seed=42)
    BA_graph = nx.read_edgelist('/home/varsha/Documents/final_results/input/WDS.edges')
    #n=20
    #w = {i: random.expovariate(5.0) for i in range(n)}
    #BA_graph = nx.geographical_threshold_graph(20, 50, weight=w)

    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_BA,
                                metrics=metrics,graph_name = 'WDS', generator = "musketeer_all_rescale",n_jobs=5)
    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_coarse, metrics=metrics,
                                graph_name='WDS', generator="musketeer_coarse_rescale", n_jobs=5)
    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_fine, metrics=metrics,
                                graph_name='WDS', generator="musketeer_fine_rescale", n_jobs=5)



def test_planar_rescale():
    num_replicas = 1
    #params_MUSKETEER = {'verbose': False, '-k':True,'edge_edit_rate': [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01],
     #                   'node_growth_rate': [0],
      #                  'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
       #                 'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_MUSKETEER = {'verbose': False, '-k': True,
                        'edge_edit_rate': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        'node_growth_rate': [0.10, 0.10, 0.10, 0.10, 0.10, 0.10,0.10, 0.10],
                        'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                        'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_coarse = {'verbose': False, '-k': True,
                        'edge_edit_rate': [ 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.10],
                        'node_growth_rate': [0.0, 0.0, 0.0, 0.0, 0.30, 0.30, 0.30, 0.30],
                        'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                        'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_fine = {'verbose': False, '-k': True,
                     'edge_edit_rate': [0.1, 0.1, 0.10,0.0, 0.0, 0.0,],
                      'edge_growth_rate': [0.30, 0.30, 0.30, 0.30,0.30, 0.0, 0.0, 0.0, 0.0],
                     'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                     'num_trial_particles': 50, 'num_insertion_trials': 50}
    metrics = graphutils.default_metrics[:]
    metrics = [met for met in metrics if
               met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree assortativity',
                                   'degree connectivity']]

    ######################
    #     3rd figure
    ######################
    print('Barabasi-Albert')
    params_BA = params_MUSKETEER.copy()
    params_BA['locality_bias_correction'] = []
    #metrics = [met for met in graphutils.default_metrics if met['name'] == 'num edges']
    metrics[-1]['optional'] = 0
    #BA_graph = nx.barabasi_albert_graph(n=20, m=3, seed=42)
    BA_graph = nx.read_edgelist('/home/varsha/Documents/sample_graph.edges')
    #n=20
    #w = {i: random.expovariate(5.0) for i in range(n)}
    #BA_graph = nx.geographical_threshold_graph(20, 50, weight=w)

    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_BA, metrics=metrics,
     graph_name = 'sample_graph_1', generator = "musketeer_all_rescale",n_jobs=1)
    #testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_coarse, metrics=metrics,
     #                          graph_name='sample_graph_1', generator="musketeer_coarse_rescale", n_jobs=1)
    #testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_fine, metrics=metrics,
     #                          graph_name='sample_graph_1', generator="musketeer_fine_rescale", n_jobs=1)

def test_planar_2():
    num_replicas = 30
    params_fine = {'verbose': False, '-k': True,
                     'edge_edit_rate': [0.1, 0.1, 0.1, 0.10,0.0, 0.0, 0.0, 0.0, 0.0],
                     'node_growth_rate': [0],
                     'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                     'num_trial_particles': 50, 'num_insertion_trials': 50}
    metrics = graphutils.default_metrics[:]
    metrics = [met for met in metrics if
               met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree assortativity',
                                   'degree connectivity']]

    BA_graph = nx.read_edgelist('/home/varsha/Documents/final_results/input/powergrid_subgraph.edges')

    testers.replica_vs_original(G=BA_graph, generator_func="RMAT", num_replicas=num_replicas, seed=1, params=params_fine,
                               metrics=metrics, graph_name='powergrid_subgraph',
                               generator="/home/varsha/Documents/final_results/RMAT/planar/powergrid/", n_jobs=10)

    BA_graph = nx.read_edgelist('/home/varsha/Documents/final_results/input/nasa_planar.edges')

    testers.replica_vs_original(G=BA_graph, generator_func = "RMAT",num_replicas=num_replicas, seed=1, params=params_fine,
                                metrics=metrics,graph_name = 'nasa_planar', generator = "/home/varsha/Documents/final_results/RMAT/planar/Boeing/",n_jobs=10)
    BA_graph = nx.read_edgelist('/home/varsha/Documents/final_results/input/sample_graph_1.edgelist')

    testers.replica_vs_original(G=BA_graph, generator_func = "RMAT",num_replicas=num_replicas,
                                seed=1, params=params_fine, metrics=metrics,graph_name = 'road_network',
                                generator = "/home/varsha/Documents/final_results/RMAT/planar/road_network_1/",n_jobs=10)
    BA_graph = nx.read_edgelist('/home/varsha/Documents/final_results/input/WDS.edges')

    testers.replica_vs_original(G=BA_graph, generator_func = "RMAT",num_replicas=num_replicas, seed=1, params=params_fine, metrics=metrics,graph_name = 'WDS', generator = "/home/varsha/Documents/final_results/RMAT/planar/WDS/",n_jobs=10)

def generate_rescaled_network():
    num_replicas = 1
    #params_MUSKETEER = {'verbose': False, '-k':True,'edge_edit_rate': [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01],
     #                   'node_growth_rate': [0],
      #                  'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
       #                 'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_all = {'verbose': False, '-k': True,
                        'edge_edit_rate': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        'node_growth_rate': [0.10, 0.10, 0.10, 0.10, 0.10, 0.10,0.10, 0.10],
                        'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                        'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_coarse = {'verbose': False, '-k': True,
                        'edge_edit_rate': [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        'node_growth_rate': [0.0, 0.0, 0.0, 0.0, 0.30, 0.30, 0.30, 0.30],
                        'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                        'num_trial_particles': 50, 'num_insertion_trials': 50}
    params_fine = {'verbose': False, '-k': True,
                     'edge_edit_rate': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                      'node_growth_rate': [0.30, 0.30, 0.30, 0.30,0.0, 0.0, 0.0, 0.0, 0.0],
                     'enforce_connected': True, 'accept_chance_edges': 1.0, 'num_pairs_to_sample': 150,
                     'num_trial_particles': 50, 'num_insertion_trials': 50}
    metrics = graphutils.default_metrics[:]
    metrics = [met for met in metrics if
               met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree assortativity',
                                   'degree connectivity']]

    BA_graph = nx.read_edgelist('/home/varsha/Documents/sample_graph.edges')
    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_all, metrics=metrics,
     graph_name = 'test_graph', generator = "musketeer_all_rescale",n_jobs=1)
    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_coarse, metrics=metrics,
     graph_name = 'test_graph', generator = "musketeer_all_coarse",n_jobs=1)
    testers.replica_vs_original(G=BA_graph, num_replicas=num_replicas, seed=1, params=params_fine, metrics=metrics,
     graph_name = 'test_graph', generator = "musketeer_fine",n_jobs=1)

if __name__ == '__main__': 
    #paper_benchmarks_and_figures()    
    #paper_findings_0()
    #paper_findings_1(seed=140371, intermediates=True, num_replicas=150)
    #paper_findings_2()
    #paper_findings_emergent()

    generate_rescaled_network()
    #test_planar_rescale()
    #paper_findings_one_metric(fname_data="precomputed_data/alternatives_seir_potterat_1__2015_03_06__18_13_03.pkl")
    #paper_findings_snapshots()
    #paper_alternatives_benchmarks()
    #paper_findings_memoriless()
    #paper_findings_epidemic()
    #paper_findings_extended()

    #poster_figure_potterat()
    #presentation_figures()
    #presentation_figures_annotation()
    #presentation_figures_epi()
    #presentation_figures_social()

    #self_dissimilarity_demo()
    #alg_alg(alg1=algorithms.generate_graph, alg2=algorithms.generate_graph, )   #this is a consistency test, since it calls the same algorithm twice


    #WS200 = nx.watts_strogatz_graph(n=200, k=4, p=0.02)
    #BA200 = nx.barabasi_albert_graph(n=200, m=10, seed=42)
    #BA100 = nx.barabasi_albert_graph(n=100, m=5, seed=42)
    #bias_benchmark(G=BA200, alg=algorithms.generate_graph, seed=None)
    #bias_benchmark(G=graphutils.load_graph('data-social/911_unwted.gml'), alg=algorithms.generate_graph, params={}, seed=None)
    #bias_benchmark(G=graphutils.load_graph('data-epi/potterat.gml'), alg=algorithms.generate_graph, params={}, seed=None)
    #pylab.show()

    #running_time_analysis()

    #editing_demo()

    #bearman_GCC = nx.read_edgelist('data-social/bearman_chains_edge_GCC.csv', delimiter=',')
    #assert bearman_GCC.number_of_nodes() == 288
    #assert bearman_GCC.number_of_edges() == 291
    #testers.replica_vs_original(G=bearman_GCC, 
    #        num_replicas=20, seed=20,
    #        params={'edge_edit_rate':[0.08, 0.07], 'node_edit_rate':[0.08, 0.07], 'verbose':False,  'locality_bias_correction':[] })
    #pylab.show()

    #testers.replica_vs_original(G=graphutils.load_graph('data-cyber-small/gr2.gml'), num_replicas=150, seed=None, params={'edge_edit_rate':[0.08, 0.07], 'node_edit_rate':[0.08, 0.07], 'verbose':False, 'enforce_connected':True, 'locality_bias_correction':[], 'edge_welfare_fraction':0})
    #testers.replica_vs_original(G=graphutils.load_graph('data-social/ftp_core.elist'), num_replicas=30, seed=None, params={'edge_edit_rate':[0.08, 0.07], 'node_edit_rate':[0.08, 0.07], 'verbose':False, 'enforce_connected':True, 'locality_bias_correction':[], 'edit_method':'alternating', 'edge_welfare_fraction':0})
    #pylab.show()


    ###a weighted network
    #testers.replica_vs_original(G=graphutils.load_graph('data-social/newman06_netscience.gml'), num_replicas=150, seed=10,
    #                    params={'edge_edit_rate':[0.05], 'verbose':False})
    
    #replica_vs_original_degree_distribution(seed=None)
   
    #movie_of_replicas()

    #dissimilarity_preservation_demo()
    #pylab.show()

    #sims_on_powergrids()
