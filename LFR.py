import subprocess
from networkit import *
import networkx as nx
import os
import random
import math
import time
import numpy as npr
timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime()) + '_%d'%npr.randint(1000)

def LFR_generator(G,params):
    G_new= graphio.readGraph('road_networks/sample_graph_2.edgelist', Format.EdgeListTabZero)
    gen = generators.LFRGenerator.fit(G_new, scale=1)
    tmp_replica = gen.generate()
    graphio.writeGraph(tmp_replica,'/home/varsha/Documents/OGDF/mycode/temp.gml', Format.GML)
    subprocess.call("/home/varsha/Documents/OGDF/example")
    os.remove("/home/varsha/Documents/OGDF/mycode/temp.gml")
    return G

def RMAT_generator(G_tmp,params):
    G = graphio.readGraph('road_networks/BA_graph_planar.edgelist', Format.EdgeListTabZero)
    scale = math.log(len(G.nodes()),2)
    print (scale)
    edgeFactor = 5
    a = round(random.uniform(0, .5),2)
    print(a)
    b= 0.5-a
    c= round(random.uniform(0, .5),2)
    print(b)
    d= 0.5-c
    gen = generators.RmatGenerator(scale, edgeFactor, a, b, c, d )
    gen.fit(G, scale=1, kronfit = False)
    G_replica = gen.generate()
    graphio.writeGraph(G_replica, '/home/varsha/Documents/OGDF/mycode/temp.gml', Format.GML)
    subprocess.call("/home/varsha/Documents/OGDF/example")
    os.remove("/home/varsha/Documents/OGDF/mycode/temp.gml")
    return G_tmp


def BTER_generator(G,params):
    files = os.listdir("/home/varsha/Desktop/feastpack_v1.1/output/road_networks/")
    for file in files:
        replica = nx.read_edgelist("/home/varsha/Desktop/feastpack_v1.1/output/road_networks/" + file)
        graphio.writeGraph(replica, '/home/varsha/Documents/OGDF/mycode/temp.gml', Format.GML)
        subprocess.call("/home/varsha/Documents/OGDF/example")



def Kronfit_generator():
    in_direc = "/home/varsha/Documents/final_results/input/"
    out_direc = "/home/varsha/Documents/final_results/Krogen/"
    files = os.listdir(in_direc)
    for file in files:
        complete_path = in_direc + file
        print(complete_path)
        for i in range(0,30):
            subprocess.call(["/home/varsha/Pictures/snap-master/examples/krongen/krongen", complete_path , out_direc+timeNow()+file])


def Kronfit_generator_2():
        out_direc = "/home/varsha/Documents/final_results/Krongen/generated/Boeing/"
        for i in range(0,30):
            print(i)
            subprocess.call(["/home/varsha/Pictures/snap-master/examples/krongen/krongen", '/home/varsha/Documents/final_results/nasa_planar.edges' , out_direc+timeNow()])



if __name__ == '__main__':
