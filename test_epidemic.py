'''
Computes the resilience of a graph to an epidemic
INTERNAL

Copyright (c) 2011-2014 by Alexander Gutfraind and Ilya Safro. 
All rights reserved.

Use and redistribution of this file is governed by the license terms in
the LICENSE file found in the project's top-level directory.


Notes:
-----

includes:
    various implementations of various measures;  invokes resMetricsX.py
    tests

no guarantee of support for directed graphs:
    the definition of neighbors might need to change

'''

import sys, os
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
#sys.path.append('./extensions')

import timeit, pdb
import numpy as np
import numpy.random as npr
import networkx as nx

def extentSIRSimL(G, params={}):
    #list-dict implementation of the contagion. much faster than a queue or a list-list
    #lavishly expensive: 5 repetitions per G.node
    sys.path.append('/home/gfriend/lib/python')
    sys.path.append('/homes/gfriend/lib/python')
    nn          = G.number_of_nodes()
    tau		    = params.get('tau', 0.5)
    repetitions = params.get('repetitions', 5)
    if 'tau' not in params:
        print('Warning: using detault tau!')

    if G.nodes() == []: return 0., 0.

    extents = np.zeros(nn*repetitions, dtype=np.double)
    for startN in G.nodes():
        for r in range(repetitions):
            extents[r+startN*r] += extentSIRSim_helper(G=G, startN=startN, params={'tau':tau}, vaccinated=[], weighted=False)
            
    return np.average(extents), np.std(extents)


def extentSIRSim_helper(G, startN, params, vaccinated, weighted=False, verbose=False):
    tau = params['tau']
    uninfectable = dict.fromkeys(vaccinated, 'vaccinated')
    uninfectable[startN] = 'infected'
   
    infected    = G.number_of_nodes()*[None]
    infected[0] = startN
    head        = 0
    tail        = 1
    if weighted:
        while head < tail:
            curNode = infected[head]
            head   += 1
            for nb in G.neighbors(curNode):
                wt = G.get_edge_data(curNode, nb)['weight']
                if not (nb in uninfectable) and tau/wt > npr.rand():
                    uninfectable[nb] = 'infected'
                    infected[tail]   = nb
                    tail            += 1
    else:
        while head < tail:
            curNode = infected[head]
            head   += 1
            for nb in G.neighbors(curNode):
                if not (nb in uninfectable) and tau > npr.rand():
                    uninfectable[nb] = 1
                    infected[tail]   = nb
                    tail            += 1
    if verbose:
        return len(uninfectable), uninfectable
    else:
        return len(uninfectable)

    #attempts at optimizations do not seem to produce results:
    #replacing infected/uninfectable with numpy arrays, replacing uninfectable with set


def extentSIRSimX(G, params):
    import epidemic_sim
    return epidemic_sim.extentSIR_adaptive(G=G, params=params)


############################

extentSIRSim = extentSIRSimL

############################


def profileCythonSim():
    params = {'tau':0.1, 'tolerance':0.05, 'minReplications':100}

    print('Python implementation:')
    G = nx.generators.erdos_renyi_graph(100,0.3)
    t = timeit.default_timer() 
    print(extentSIRSimL(G, params=params))
    print(timeit.default_timer()-t)

    print('Cython implementation:')
    t = timeit.default_timer() 
    try:
      print(extentSIRSimX(G, params=params))
    except Error as e:
      print(e)
    print(timeit.default_timer()-t)



def testCython():
    import pylab

    meansCyt = []
    sdevsCyt = []
    meansNum = []
    sdevsNum = []
    params = {'minReplications':5, 'tolerance':0.1, 'tau':0.1}
    for graphNum in range(10):
        G = nx.generators.erdos_renyi_graph(100,0.2)
        for trial in range(10):
            m, s  = extentSIRSimX(G, params)
            meansCyt.append(m)
            sdevsCyt.append(s)
            m, s  = extentSIRSim(G, params)
            meansNum.append(m)
            sdevsNum.append(s)
    pylab.show()
    figMC = pylab.figure()
    pylab.hold(True)
    bins = list(range(0,101,10))
    countsMeanCyt = pylab.hist(meansCyt, bins=bins, histtype='step', label='Cyt')[0]
    countsMeanNum = pylab.hist(meansNum, bins=bins, histtype='step', label='Num')[0]
    pylab.title('means') 

    print()
    print('Means')
    meansDiff = sum(np.abs(countsMeanCyt-countsMeanNum))
    relError = 0.  #aside: rel error increases when tau decreases -> smaller extents
    for i, mn in enumerate(meansCyt):
        if mn > 0.:
            relError += (mn-meansNum[i])/mn
    print('Absolute difference (always>0): %f'%meansDiff)
    print('Avg relative error: %f'%(relError/len(meansCyt)))
    #print 'Net      difference: %f'%sum(      (countsMeanCyt-countsMeanNum))
    if (meansDiff / float(len(bins)*(graphNum+1))) < 1 and relError < 0.2:
        print('Test passed!')
        print('Visually inspect: the histograms should broadly overlap')
    else:
        print('Test failed: more than one misplaced count per graph per bin!')


    print()
    print('SDevs')
    pylab.show()
    figMC = pylab.figure()
    pylab.hold(True)
    bins = list(range(0,101,10))
    countsStdCyt = pylab.hist(sdevsCyt, bins=bins, histtype='step', label='Cyt')[0]
    countsStdNum = pylab.hist(sdevsNum, bins=bins, histtype='step', label='Num')[0]
    pylab.title('sdevs') 

    stdevsDiff = sum(np.abs(countsStdCyt-countsStdNum))
    relError = 0.
    for i, sd in enumerate(sdevsCyt):
        if sd > 0.:
            relError += (sd-sdevsNum[i])/mn
    print('Absolute difference: %f'%stdevsDiff)
    print('Avg of relative error: %f'%(relError/len(sdevsCyt)))


def testCythonExtentLong():
    import pylab
    import epidemic_sim
    params = {'minReplications':1, 'tau':0.4, 'tie_stability':0.5, 'latent_interval':3, 'infectious_interval':2}
    G = nx.generators.erdos_renyi_graph(500,0.02)
    
    figMC = pylab.figure()
    pylab.hold(True)
    
    for trial in range(10):
        mean_extents = epidemic_sim.extentSIR_long(G=G, params=params)
        print(sum(mean_extents))
        if sum(mean_extents) > 0:
            pylab.plot(mean_extents[:1+max(mean_extents.nonzero()[0])])

    pylab.hold(False)
    pylab.xlabel('time (day)')
    pylab.ylabel('new cases')
    pylab.show()

if __name__ == '__main__':
    print('Running Tests...\n')
    #profileCythonSim()
    #testCython()

    testCythonExtentLong()
