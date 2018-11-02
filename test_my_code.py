import planarity
import networkx as nx
import Algebraic_Rank
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy.io import mminfo,mmread
import random
from random import choice

G = nx.read_edgelist("mesh33.edges")
for edge in G.edges():
       G.add_edge(edge[0], edge[1], weight=1)
nb=0
added_edges=[]
while nb < 5:
    first_node = choice(G.nodes())  # pick a random node
    possible_nodes = set(G.nodes())
    neighbours = G.neighbors(first_node) + [first_node]
    possible_nodes.difference_update(neighbours)  # remove the first node and all its neighbours from the candidates
    second_node = choice(list(possible_nodes))
    G.add_edge(first_node, second_node, weight=1)
    added_edges.append((first_node,second_node))
    added_edges.append((second_node,first_node))
    print(first_node, second_node)
    nb += 1
pos = graphviz_layout(G, prog='sfdp')
nx.draw(G, pos, with_labels=False, node_size=1)
plt.show()


itr = 7
E_final = {}

for edge in G.edges():
    E_final[edge] = 0
while itr > 0:
    X = Algebraic_Rank.get_algebraic_rank(G, 0.5, 100)
    E = {}
    for edge in G.edges():
        E[edge] = abs(X[edge[0]] - X[edge[1]])
        if(E_final[edge]<E[edge]):
            E_final[edge]=E[edge]
    itr -= 1
itr=0
E_sorted_reverse=[]
E_sorted = OrderedDict(sorted(E_final.items(), key=lambda x: x[1], reverse=True))
for edge in E_sorted:
        E_sorted_reverse.append((edge[1],edge[0]))
        if edge in added_edges:
            print(edge , itr)
        itr+=1

nb_of_edges_removed=0

print(nb_of_edges_removed,planarity.is_planar(G.edges()))
while not planarity.is_planar(G.edges()):
        K_edges = planarity.kuratowski_edges(G.edges())
        K_edges_reverse = []
        for edge in K_edges:
            K_edges_reverse.append((edge[1],edge[0]))
        for edge in E_sorted:
            if edge not in K_edges and edge not in K_edges_reverse:
               continue
            else:
               G.remove_edge(edge[0], edge[1])
               nb_of_edges_removed += 1
               break

   # print(itr)
pos = graphviz_layout(G, prog='sfdp')
nx.draw(G, pos, with_labels=False, node_size=1)
plt.show()
print(nb_of_edges_removed,planarity.is_planar(G.edges()))
