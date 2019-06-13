import networkx as nx
import planarity


def make_planar(G):
    if planarity.is_planar(G.edges()):
        return G
    all_kuratowski_edges = []
    next_G = G.copy()

    while planarity.is_planar(G.edges()) == False:
        edge_list = planarity.kuratowski_edges(G)
        print(len(edge_list))
        G.remove_edges_from(edge_list)
        for edge in edge_list:
            all_kuratowski_edges.append(edge)

    itr=0
    replaced = []
    for edge in all_kuratowski_edges:
        neighbors = nx.neighbors(G,edge[1])
        if itr < len(neighbors) and edge not in replaced:
            G.add_edge(edge[0],neighbors[itr])
        elif itr >= len(neighbors) and edge not in replaced:
            global neighbors
            neighbors = nx.neighbors(G,edge[0])
            itr=0
        itr = itr+1


    print (planarity.is_planar(G),nx.number_of_edges(G))
    print (planarity.is_planar(next_G),nx.number_of_edges(next_G))


