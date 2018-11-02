from networkit import *
import networkx as nx
from random import *
import os, subprocess
from datetime import datetime
import time
import math
import random

def LFR_generator(G,path):
    gen = generators.LFRGenerator.fit(G, scale=1)
    G_replica = gen.generate()
    graphio.writeGraph(G_replica, path, Format.GML)
    return G_replica

def BTER_generator(G):
    gen = generators.BTERReplicator.fit(G, scale = 1)
    gen.setPaths(feastpackPath='/home/varsha/Desktop/feastpack_v1.1')
    G_replica = gen.generate()
    return G_replica

def RMAT_generator(G,path):
    scale = math.log(len(G.nodes()), 2)
    print(scale)
    edgeFactor = 5
    a = round(random.uniform(0, .5), 2)
    print(a)
    b = 0.5 - a
    c = round(random.uniform(0, .5), 2)
    print(b)
    d = 0.5 - c
    gen = generators.RmatGenerator(scale, edgeFactor, a, b, c, d)
    gen.fit(G, scale=1, kronfit=False)
    G_replica = gen.generate()
    graphio.writeGraph(G_replica, path, Format.GML)
    return G_replica

def test_compare():
    mesh_G = graphio.readGraph('road_networks/mesh33.edgelist', Format.EdgeListTabZero)
    gen = generators.LFRGenerator.fit(mesh_G, scale=1)
    mesh_G_replica = gen.generate()
    graphio.writeGraph(mesh_G_replica,'road_networks/mesh33_replica.gml', Format.GML)

def other_generators():
    test_planar_1(params_coarse_changes, mesh_G)
    test_planar_1(params_fine_changes, mesh_G)
    test_planar_1(params_all_level_changes, mesh_G)

    RN_1 = nx.read_edgelist('road_networks/sample_graph_1.edges')
    test_planar_1(params_coarse_changes, RN_1)
    test_planar_1(params_fine_changes, RN_1)
    test_planar_1(params_all_level_changes, RN_1)

    RN_2 = nx.read_edgelist('road_networks/sample_graph_2.edges')
    test_planar_1(params_coarse_changes, RN_2)
    test_planar_1(params_fine_changes, RN_2)
    test_planar_1(params_all_level_changes, RN_2)

    WDS = nx.read_edgelist('road_networks/WDS.edges')
    test_planar_1(params_coarse_changes, WDS)
    test_planar_1(params_fine_changes, WDS)
    test_planar_1(params_all_level_changes, WDS)

    BA_G = nx.read_gml('road_networks/BA_subgraph.gml')
    test_planar_1(params_coarse_changes, BA_G)
    test_planar_1(params_fine_changes, BA_G)
    test_planar_1(params_all_level_changes, BA_G)

    ER_G = nx.read_gml('road_networks/ER_subgraph.gml')
    test_planar_1(params_coarse_changes, ER_G)
    test_planar_1(params_fine_changes, ER_G)
    test_planar_1(params_all_level_changes, ER_G)

    Cmp_G = nx.read_gml('road_networks/Cmp_subgraph.gml')
    test_planar_1(params_coarse_changes, Cmp_G)
    test_planar_1(params_fine_changes, Cmp_G)
    test_planar_1(params_all_level_changes, Cmp_G)

    PG_G = nx.read_gml('road_networks/powergrid_subgraph.gml')
    test_planar_1(params_coarse_changes, PG_G)
    test_planar_1(params_fine_changes, PG_G)
    test_planar_1(params_all_level_changes, PG_G)

    Random_G = nx.read_gml('road_networks/Random_planar_subgraph.gml')
    test_planar_1(params_coarse_changes, Random_G)
    test_planar_1(params_fine_changes, Random_G)
    test_planar_1(params_all_level_changes, Random_G)

    Rnd_G = nx.read_gml('road_networks/Rnd_reg_graph.gml')
    test_planar_1(params_coarse_changes, Rnd_G)
    test_planar_1(params_fine_changes, Rnd_G)
    test_planar_1(params_all_level_changes, Rnd_G)


path_1 = '/home/varsha/Documents/final_results/LFR/generated/powergrid/'
path_2 = '/home/varsha/Documents/final_results/RMAT/generated/powergrid/'
G = graphio.readGraph('/home/varsha/Documents/final_results/input/powergrid_subgraph.edges', Format.EdgeListSpaceOne)
for i in range(0,30):
    time = str(datetime.now().time()) + ".gml"
    LFR_generator(G,path_1+ time)
    subprocess.call(['/home/varsha/Documents/OGDF/planarize',path_1+time])
    RMAT_generator(G, path_2 + time)
    subprocess.call(['/home/varsha/Documents/OGDF/planarize', path_2 + time])
