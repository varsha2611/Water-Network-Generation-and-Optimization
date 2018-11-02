import planarity
import networkx as nx
import Algebraic_Rank
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy.io import mminfo,mmread
import random
from random import choice
def nodes_in_K_graph(K_edges):
    nodes = {}
    for edge in K_edges:
        nodes[edge[0]] = 0
        nodes[edge[1]] = 0
    for edge in K_edges:
        nodes[edge[0]]+=1
        nodes[edge[1]]+=1
    nodes_sorted = OrderedDict(sorted(nodes.items(), key=lambda x: x[1], reverse=True))
    main_nodes=[]
    for node in nodes_sorted:
        if nodes_sorted[node]>2:
            main_nodes.append(node)
        else:
            break
    return main_nodes

G = nx.read_edgelist("mesh33.edges")
for edge in G.edges():
       G.add_edge(edge[0], edge[1], weight=1)
nb=0
added_edges=[]
added_edges_node =[]
while nb < 5:
    first_node = choice(G.nodes())  # pick a random node
    possible_nodes = set(G.nodes())
    neighbours = G.neighbors(first_node) + [first_node]
    possible_nodes.difference_update(neighbours)  # remove the first node and all its neighbours from the candidates
    second_node = choice(list(possible_nodes))
    G.add_edge(first_node, second_node, weight=1)
    added_edges.append((first_node,second_node))
    added_edges.append((second_node,first_node))
    added_edges_node.append(first_node)
    added_edges_node.append(second_node)
    print(first_node, second_node)
    nb += 1
pos = graphviz_layout(G, prog='sfdp')
nx.draw(G, pos, with_labels=False, node_size=1)
plt.show()
print(len(G.nodes()),planarity.is_planar(G.edges()))
itr =1
while not planarity.is_planar(G.edges()):
        K_edges = planarity.kuratowski_edges(G.edges())
        nodes = nodes_in_K_graph(K_edges)
        i=0
        for node in nodes:
            for edge in G.edges(node):
                if edge in K_edges or (edge[1],edge[0]) in K_edges:
                    G.remove_edge(edge[0], edge[1])
                    removed_edge.append((edge[0],edge[1]))

        itr+=1

print(itr)

nb_of_edges_removed = 0
for edge in removed_edge:
    G.add_edge(edge[0],edge[1])
    if planarity.is_planar(G):
        continue
    else:
        G.remove_edge(edge[0],edge[1])

pos = graphviz_layout(G, prog='sfdp')
nx.draw(G, pos, with_labels=False, node_size=1)
plt.show()
print(len(G.nodes()),planarity.is_planar(G.edges()))
