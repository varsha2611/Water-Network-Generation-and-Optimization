import networkx as nx
import collections
from collections import deque
from random import choice

import networkx as nx
from collections import defaultdict, deque

def bfs_edges(G, n, source, reverse=False):
    if reverse and isinstance(G, nx.DiGraph):
        neighbors = G.predecessors_iter
    else:
        neighbors = G.neighbors_iter
    visited=set([source])
    queue = deque([(source, neighbors(source))])
    while queue and len(visited)<n:
        parent, children = queue[0]
        try:
            child = next(children)
            if child not in visited:
                yield parent, child
                visited.add(child)
                queue.append((child, neighbors(child)))
        except StopIteration:
            queue.popleft()

def bfs_tree_custom(G,n):
    source = choice(G.nodes())
    print "source",source
    T = nx.Graph()
    T.add_node(source)
    T.add_edges_from(bfs_edges(G,n,source,reverse=False))
    return T

def bfs_predecessors(G, source):
    return dict((t,s) for s,t in bfs_edges(G,source))

def bfs_successors(G, source):
    d=defaultdict(list)
    for s,t in bfs_edges(G,source):
        d[s].append(t)
    return dict(d)



