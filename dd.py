# -*- coding: utf-8 -*-
'''
Multiscale Entropic Network Generator (MUSKETEER)

Copyright (c) 2011-2014 by Alexander Gutfraind and Ilya Safro. 
All rights reserved.

Use and redistribution of this file is governed by the license terms in
the LICENSE file found in the project's top-level directory.

analysis of deferential detachment

'''

import os
import time
import numpy as np
import numpy.random as npr
import random, sys
import networkx as nx
import matplotlib
#matplotlib.use('PS')
#matplotlib.use('PDF')
#import matplotlib.pylab as pylab
import pylab
import pdb
import pickle
import scipy.stats

#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#matplotlib.rc('font',**{'family':'serif','serif':['Palatino']})
#rc('font',**{'family':'serif','serif':['Times New Roman']})
matplotlib.rc('text', usetex=True)

np.seterr(all='raise')

timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())

def compare_randomization():
    G = nx.erdos_renyi_graph(100, 0.1)
    #G = nx.barabasi_albert_graph(100, 10)

    base_seq = list(nx.degree(G).values())
    base_seq.sort()
    
    num_edits = 10000
    dd_seq   = run_dd2(base_seq[:], num_edits)
    dd_seq.sort()
    edit_seq = run_random_edits(base_seq[:], num_edits)
    edit_seq.sort()

    pylab.figure()
    pylab.hold(True)
    pylab.plot(base_seq, 'r-o', label='original')
    pylab.plot(dd_seq,   'b-D', label='dd')
    pylab.plot(edit_seq, 'y-*', label='random')

    pylab.legend()
    pylab.savefig('output/randomization_sequences_'+timeNow())
    pylab.show()

def run_dd1(base_seq, num_edits):
    new_seq = base_seq[:]
    num_nodes = len(base_seq)
    for i in range(num_edits):
        loser_node = npr.randint(num_nodes)
        d = new_seq[loser_node]
        if d > 0:
            new_seq[loser_node] -= npr.binomial(d, 1./d)
            new_seq[npr.randint(num_nodes)] += 1

    return new_seq

def run_dd2(base_seq, num_edits):
    new_seq = base_seq[:]
    num_nodes = len(base_seq)
    for i in range(num_edits):
        loser_node = npr.randint(num_nodes)
        d = new_seq[loser_node]
        if d > 0 and npr.rand() < 1./d:
            new_seq[loser_node] -= 1
            new_seq[npr.randint(num_nodes)] += 1

    return new_seq

def run_random_edits(base_seq, num_edits):
    new_seq = base_seq[:]
    num_nodes = len(base_seq)
    for i in range(num_edits):
        loser_node = npr.randint(num_nodes)
        if new_seq[loser_node] > 0:
            new_seq[loser_node] -= 1
            new_seq[npr.randint(num_nodes)] += 1

    return new_seq

if __name__ == '__main__': 
    compare_randomization()
