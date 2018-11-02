'''
Computes the extent of an epidemic (implementation in Cython)

internal

Notes:
-----

'''

import sys
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
#sys.path.append('./extensions')

import numpy as np
cimport numpy as np
import numpy.random as npr

import networkx as nx 

cdef extern from "stdlib.h":
    #TODO: random,srandom fails to compile under mingw.  does new code work in linux?
    #long c_libc_random "random"()
    #void c_libc_srandom "srandom"(unsigned int seed)
    long c_libc_random "rand"()
    void c_libc_srandom "srand"(unsigned int seed)
    void free(void* ptr)
    void* malloc(size_t size)
    void* realloc(void* ptr, size_t size)
    
        
#from stdlib cimport *

DEF MAX_SEED = 32767

cpdef extentSIR_long(G, dict params={}):
#reports the temporal structure (incidence over time) and final size distribution of the final number of infectees
    cdef int nn              = G.number_of_nodes()
    cdef float tau	         = params.get('tau', 0.5)
    cdef float tie_stability = params.get('tie_stability', 0.9)
    cdef int latent_interval       = params.get('latent_interval', 3)
    cdef int infectious_interval   = params.get('infectious_interval', 7)
    cdef int max_time        = params.get('max_time', nn*(infectious_interval+latent_interval))

    cdef int minReplications = params.get('minReplications', max(nn/10, 40))

    if nn==0: return 0., 0.  #this is a problem since we assume initial infection site
    elif nn==1: return 1., 0. 

    G_int = nx.convert_node_labels_to_integers(G)

    c_libc_srandom(npr.randint(low=MAX_SEED)) #could potentially pass the seed and avoid importing npr
    
    cdef int **lol = <int **>malloc(nn * sizeof(int*)) #graph as list-of-lists
    cdef int node,
    cdef list nbs
    cdef int numNbs
    cdef int j
    for node in range(nn): #range is apparently better than xrange
        nbs = G_int.neighbors(node)
        #print nbs
        numNbs = len(nbs)
        lol[node] = <int *>malloc((1+numNbs) * sizeof(int))
        lol[node][0] = numNbs #0 stores the number of neighbors
        for j from 1<=j<=numNbs:
            lol[node][j] = nbs[j-1]

    cdef int trial   = 0
    cdef int startNode
    #print minReplications
    cdef np.ndarray[np.int_t, ndim=1] mean_infected = np.zeros(max_time, dtype=int)
    cdef np.ndarray[np.int_t, ndim=1] sizes         = np.zeros(nn+1, dtype=int)
    cdef np.ndarray[np.int_t, ndim=1] new_infected
    cdef int total_infected
    while trial < minReplications:
        #print trial
        startNode = c_libc_random() % nn
        #print
        #print startNode
        new_infected = np.zeros(max_time, dtype=int)
        #print len(new_infected)
        runSingleNodeSEIR(lol, nn, startNode, tau, tie_stability, infectious_interval, latent_interval, new_infected)
        mean_infected += new_infected
        total_infected = new_infected.sum()
        sizes[total_infected] += 1
        trial    += 1

    #print 'Num trials: %d'%trial 

    #print
    #print sizes
    #print mean_infected
    for node in range(nn):
        free(lol[node])
    free(lol)
    return (mean_infected / float(minReplications), sizes / float(minReplications))

cdef void runSingleNodeSEIR(int** lol, int nn, int startN, float tau, float tie_stability, int infectious_interval, int latent_interval, np.ndarray[np.int_t, ndim=1] new_infected):
    #wishlist: add variance in latency, infectious period
    #cdef int intThreshold = <int>(MAX_SEED*tau)
    cdef set uninfectable = set()
    cdef dict infected    = {startN:infectious_interval}  #node->remaining days
    cdef dict latent      = {}                            #node->remaining days

    #cdef np.ndarray[np.int_t, ndim=1] new_infected  = np.zeros(nn*LATENT_PERIOD, dtype=int) #i.e. become exposed even if not infectious themselves
    #cdef np.ndarray[np.int_t, ndim=1] new_infected = np.zeros(nn*(infectious_interval+latent_interval), dtype=int)
    cdef int t = 0
 
    cdef int nb
    cdef int nbCounter
    cdef int numNbs
    #print '  ' + str(len(new_infected))
    while len(infected) + len(latent) > 0 and t < len(new_infected):
        #print 't: %d'%t
        #print 'infected: ' + str(infected)
        #print 'latent:  '  + str(latent)
        for node in infected.keys():
            if infected[node] == 0:
                infected.pop(node)
                continue
            else:
                infected[node] -= 1
            numNbs = lol[node][0]  #stored in this data structure
            for nbCounter from 1<=nbCounter<=numNbs:
                nb = lol[node][nbCounter]
                #if (nb not in uninfectable) and intThreshold > tie_stability*(c_libc_random()%MAX_SEED):
                if (nb not in uninfectable) and tau*tie_stability > npr.rand():
                    #print 'new inf: %d -> %d'%(node,nb)
                    uninfectable.add(nb)
                    latent[nb]   = latent_interval + npr.randint(-1,2)
                    new_infected[t] += 1
        for node in latent.keys():
            if latent[node] == 1:
                latent.pop(node)
                infected[node] = infectious_interval + npr.randint(-1,2)
            else:
                latent[node] -= 1

        t += 1
        #print 'time: %d'%t
    return


