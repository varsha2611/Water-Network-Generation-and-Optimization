# -*- coding: utf-8 -*-
'''
Multiscale Entropic Network Generator (MUSKETEER)

Copyright (c) 2011-2014 by Alexander Gutfraind and Ilya Safro. 
All rights reserved.

Use and redistribution of this file is governed by the license terms in
the LICENSE file found in the project's top-level directory.

Benchmarks of Quality and performance - internal 

'''

import os, subprocess
import time
import numpy as np
import numpy.random as npr
import random, sys
import networkx as nx
import matplotlib
#matplotlib.use('PS')
matplotlib.use('PDF')
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

def paper_benchmark_graphs(seed=8):
    if seed==None:
        seed = npr.randint(1E6)
    print('rand seed: %d'%seed)
    npr.seed(seed)
    random.seed(seed)

    num_replicas = 30

    params_default = {'verbose':False, 
                      'edge_edit_rate':[0.01, 0.01, 0.01, 0.01, 0.01], 'node_edit_rate':[0.01, 0.01, 0.01, 0.01, 0.01], 'node_growth_rate':[0],  
                      'enforce_connected':True, 'locality_bias_correction':[0,],
                      #'num_v_cycles':1,
                      #'algorithm': (algorithms.musketeer_on_subgraphs, algorithms.musketeer_iterated_cycle,)
                     }
    params_for_big = params_default.copy()
    metrics_default = graphutils.default_metrics[:]
    #metrics_default  = filter(lambda met: met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree connectivity', 'degree assortativity',  'average shortest path', 'mean ecc', 'powerlaw exp', ], 
    #                metrics_default)
    metrics_default = [k for k in graphutils.default_metrics if k['name'] in ['num nodes', 'num edges', 'clustering', 'total deg*deg', 'harmonic mean path']]
    #metrics_default  = filter(lambda met: met['name'] not in ['avg flow closeness', 'avg eigvec centrality', 'degree connectivity', 'degree assortativity',  'average shortest path', 'mean ecc', 'powerlaw exp', ], metrics_default)
    problems =  [{'graph_path':'data-benchmark/911_unwted.edges'},
                 {'graph_path':'data-benchmark/BA5000_10.edges'},
                 {'graph_path':'data-benchmark/ER5000_0d001.edges'},
                 {'graph_path':'data-benchmark/Fly_protein_interactions_IntAct.elist'},
                 {'graph_path':'data-benchmark/Yeast_protein_interactions_IntAct.elist'},
                 {'graph_path':'data-benchmark/as-caida20071105.elist'},
                 {'graph_path':'data-benchmark/darwinBookInter_st.elist'},
                 {'graph_path':'data-benchmark/facebook-links.elist'},
                 {'graph_path':'data-benchmark/openflights.elist'},
                 {'graph_path':'data-benchmark/p2p-Gnutella24.elist'},
                 {'graph_path':'data-benchmark/potterat_Hiv250.elist'},
                 {'graph_path':'data-benchmark/salganic-proj90-epi-gcc.welist'},
                 {'graph_path':'data-benchmark/twittercrawl.elist'},
                 {'graph_path':'data-benchmark/ukroads.adjlistImp'},
                 {'graph_path':'data-benchmark/watts_strogatz98_power.elist'},
                ]

    return_data = {}
    for problem_info in problems:
        graph_path      = problem_info['graph_path']
        graph_params    = problem_info.get('graph_params', params_default)
        graph_metrics   = problem_info.get('graph_metrics', metrics_default)
        graph_num_replicas = problem_info.get('num_replicas', num_replicas)
        #graph_num_replicas = 2
        print('Starting:')
        print(graph_path)
        base_graph      = graphutils.load_graph(path=graph_path, params={'verbose':False})
        print('Size nn: %d, ne: %d'%(nx.number_of_nodes(base_graph), nx.number_of_edges(base_graph)))
        base_graph = nx.Graph(base_graph)

        base_graph.name = graph_path.split('/')[1]

        gpath     = 'output/'+base_graph.name+'_'+timeNow()+'.dot'
        gpath_fig = gpath[:-3]+'eps'
        graphutils.write_graph(G=base_graph, path=gpath)
        try:
            ret=testers.replica_vs_original(G=base_graph, num_replicas=graph_num_replicas, seed=seed, params=graph_params, metrics=graph_metrics, title_infix='', n_jobs=-1)
            return_data[graph_path] = ret
        except Exception as inst:
            print('Warning: could not generate or evaluate '+graph_path + ': '+str(inst))
            exc_traceback = sys.exc_info()[2]
            print(str(inst) + "\n" + str(traceback.format_tb(exc_traceback)).replace('\\n', '\n'))

    graphutils.safe_pickle(path='output/sweep_report_'+timeNow()+'.pkl', data=return_data)

    print('Summary')
    print('Graph : MeanMeanError,MeanMeanErrorStd')
    mm1,mm2=[],[]
    for graph_path, data in list(return_data.items()):
        v1,v2 = data[0]['mean_mean_error'], data[0]['mean_mean_errorstd']
        print('%s : %.3f, %.3f'%(graph_path, v1, v2))
        mm1.append(v1)
        mm2.append(v2)
    print('OVERALL MEAN : %.3f, %.3f'%(np.average(mm1),np.average(mm2)))

    print()
    print('MeanMeanError = the relative error of each replica to the original is computed for each metric, and averaged over all replicas and then all metrics')
    print('MeanMeanErrorStd  = same input to above, but the mean error per metric is divided by the standard deviation of errors for that metric, and then averaged over metrics')
    

if __name__ == '__main__': 
    paper_benchmark_graphs()    


