import planarity
import networkx as nx
import algebraic_distance
import Algebraic_Rank
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy.io import mminfo,mmread
from random import choice
import time


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
    #print('no of main nodes',len(main_nodes))
    return main_nodes

def make_planar(G):
    if planarity.is_planar(G.edges()):
        print ("Generated Graph is planar")
        return G
    removed_edge =[]
    itr = 0
    G_copy = G.copy()
    for edge in G_copy.edges():
        G_copy.add_edge(edge[0], edge[1], weight=1)
    print(len(G.edges()), planarity.is_planar(G.edges()))
    while not planarity.is_planar(G):
        K_edges = planarity.kuratowski_edges(G)
        nodes = nodes_in_K_graph(K_edges)
        for edge in K_edges:
            if edge[0] in nodes or edge[1] in nodes:
                G.remove_edge(edge[0], edge[1])
                removed_edge.append((edge[0], edge[1]))
        itr += 1
    print("iterations",itr)
    print("removed edges", len(removed_edge))
    itr = 7
    E_final={}
    for edge in removed_edge:
        E_final[edge] = 0
    while itr > 0:
        X = Algebraic_Rank.get_algebraic_rank(G_copy, 0.5, 100)
        E = {}
        for edge in removed_edge:
            E[edge] = abs(X[edge[0]] - X[edge[1]])
            if (E_final[edge] < E[edge]):
                E_final[edge] = E[edge]
        itr -= 1

    E_sorted = OrderedDict(sorted(E_final.items(), key=lambda x: x[1], reverse=False))
    for edge in E_sorted:
        G.add_edge(edge[0], edge[1])
        if planarity.is_planar(G):
            del E_sorted[edge]
            continue
        else:
            G.remove_edge(edge[0], edge[1])

    start = time.time()

    print("before rewiring", len(removed_edge))
    for edge in E_sorted:
        neigbors = G.neighbors(edge[1])
        for neighbor in neigbors:
            if not G.has_edge(edge[0], neighbor) and neighbor != edge[0]:
                G.add_edge(edge[0],neighbor)
                if planarity.is_planar(G):
                    del E_sorted[edge]
                    break
                else:
                    G.remove_edge(edge[0], neighbor)
                    continue
            else:
                continue

    for edge in E_sorted:
        neigbors = G.neighbors(edge[0])
        for neighbor in neigbors:
            if not G.has_edge(edge[1], neighbor) and neighbor != edge[1]:
                G.add_edge(edge[1], neighbor)
                if planarity.is_planar(G):
                    del E_sorted[edge]
                    break
                else:
                    G.remove_edge(edge[1], neighbor)
                    continue
            else:
                continue
    end = time.time()
    # run your code
    elapsed = end - start
    print("after first level rewiring", len(E_sorted))
    print('time for first level rewiring: ', elapsed / 10, 'seconds')

    start = time.time()
    for edge in E_sorted:
        edge_rewired = 0
        neigbors = G.neighbors(edge[0])
        for neighbor in neigbors:
            neighbors_lvl_2 = G.neighbors(neighbor)
            for neighbor_lvl_2 in neighbors_lvl_2:
                if not G.has_edge(edge[1], neighbor_lvl_2) and neighbor_lvl_2 != edge[1]:
                     G.add_edge(edge[1], neighbor_lvl_2)
                     if planarity.is_planar(G):
                         del E_sorted[edge]
                         edge_rewired+=1
                         break
                     else:
                         G.remove_edge(edge[1], neighbor_lvl_2)
                         continue
                else:
                    continue

            if edge_rewired > 0:
                break

    print("after first second level rewiring", len(E_sorted))

    for edge in E_sorted:
        neigbors = G.neighbors(edge[1])
        edge_rewired = 0
        for neighbor in neigbors:
            neighbors_lvl_2 = G.neighbors(neighbor)
            for neighbor_lvl_2 in neighbors_lvl_2:
                if not G.has_edge(edge[0], neighbor_lvl_2) and neighbor_lvl_2 != edge[0]:
                     G.add_edge(edge[0], neighbor_lvl_2)
                     if planarity.is_planar(G):
                        del E_sorted[edge]
                        edge_rewired += 1
                        break
                     else:
                        G.remove_edge(edge[0], neighbor_lvl_2)
                        continue
                else:
                    continue

            if edge_rewired > 0:
                break

    end = time.time()
    # run your code
    elapsed = end - start

    print('time for second level rewiring: ', elapsed/10,'seconds')
    print("after second level rewiring", len(E_sorted))
    print(len(G.edges()), planarity.is_planar(G.edges()))
    nx.write_adjlist(G, "Generated_planar_graph.edges")
    pos = graphviz_layout(G, prog='sfdp')
    nx.draw(G, pos, with_labels=False, node_size=1)
    plt.show()
    return G
