'''
Multiscale Entropic Network Generator (MUSKETEER)

Copyright (c) 2011-2014 by Alexander Gutfraind and Ilya Safro. 
All rights reserved.

Use and redistribution of this file is governed by the license terms in
the LICENSE file found in the project's top-level directory.


Main user program

- processes user input and reads data
- calls algorithms
- returns replicas

'''

import os
import time
import numpy as np
import numpy.random as npr
import random, sys
import networkx as nx
import matplotlib
#matplotlib.use('PDF')
#import matplotlib.pylab as pylab
import pylab
import getopt
import re
import pdb
import pickle
import algorithms
import graphutils
import simpletesters

np.seterr(all='raise')
timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())

version = 'Beta 1'

G=nx.read_edgelist(sys.argv[1])
H=nx.read_edgelist(sys.argv[2])

X = nx.disjoint_union(G, H)

C1=nx.connected_component_subgraphs(X)[0]
C2=nx.connected_component_subgraphs(X)[1]

n1 = random.choice(C1.nodes())
n2 = random.choice(C2.nodes())
X.add_edge(n1, n2)

r = random.uniform(0, 1)
if r < 0.5:
    n1 = random.choice(X.nodes())
    n2 = random.choice(X.nodes())
    X.add_edge(n1, n2)

nx.write_edgelist(X, "outgraph.elist")
graphutils.write_graph(X, "outgraph.dot")
