from epynet import *
from epanettools.epanettools import EPANetSimulation,Node, Link, Pattern
import random
def SetValues(inputFile,id = "",ep = True):
    d = Network(inputFile)
    curves = len(d.curves)
    max_elevation = 0

    if(ep):
        d.add_curve('C1', [(1000, 180)])
        d.add_curve('C2', [(2000, 200)])
        d.add_curve('C3', [(3000, 220)])
        d.add_curve('C4', [(4000, 240)])
        d.add_curve('C5', [(5000, 260)])
        d.add_curve('C6', [(1000, 200)])
        d.add_curve('C7', [(2000, 220)])
        d.add_curve('C8', [(3000, 240)])
        d.add_curve('C9', [(4000, 260)])
        d.add_curve('C10', [(5000, 280)])

    MinHead =[40]
    nodes = d.nodes
    nnodes = len(nodes)
    hStar = [0]*nnodes
    nbOfPipes = 0
    nbOfPumps = 0
    nbOfTanks = 0
    links = d.links
    curves = len(d.curves)
    for link in links:
        if link.link_type == 'pipe':
            nbOfPipes+=1
        elif link.link_type == 'pump':
            nbOfPumps+=1
    nbofJunctions = 0
    for node in nodes:
        if node.node_type == 'Junction':
            hStar[nbofJunctions] = random.choice(MinHead)
            max_elevation = max(max_elevation,node.elevation)
            nbofJunctions += 1
        elif node.node_type == 'Tank':
            nbOfTanks+=1
    #d.save_inputfile('temp.inp')

    if ep:
        d.save_inputfile('temp'+id+'.inp')
        et = EPANetSimulation('temp'+id+'.inp')
        Conn, NoConn = GetConnectionDetails(et)
        return et,hStar,curves, curves,nbOfPipes,nbOfPumps,nbOfTanks,Conn,NoConn,max_elevation
    else:
        return nbOfPipes,nbOfPumps,nbOfTanks,max_elevation


def GetConnectionDetails(d):
    ret, nnodes = d.ENgetcount(d.EN_NODECOUNT)
    ret, llinks = d.ENgetcount(d.EN_LINKCOUNT)
    Conn = []
    NoConn = []
    h_degree = 0
    if (len(Conn) == 0):
        for n in range(0, nnodes + 1):
            c = []
            Conn.append(c)

        # generate mapping used later for resilience calculation
    for i in range(0, (llinks + 1)):
        if (len(Conn) == (nnodes + 1)):
            nodes = d.ENgetlinknodes(i + 1)
            if (nodes[0] == 0):
                Conn[nodes[1]].append(i + 1)
                Conn[nodes[2]].append(i + 1)
                h_degree = max(len(Conn[nodes[0]]), len(Conn[nodes[1]]), h_degree)


    for idx in range(0, nnodes + 1):
        NoConn.append(len(Conn[idx]))
        while (len(Conn[idx]) < h_degree):
            Conn[idx].append(0)
    del Conn[0]
    del NoConn[0]

    return Conn,NoConn


def GetConnectionData(inputFile):
    d = EPANetSimulation(inputFile)
    ret, nnodes = d.ENgetcount(d.EN_NODECOUNT)
    ret, llinks = d.ENgetcount(d.EN_LINKCOUNT)
    Conn = []
    NoConn = []
    h_degree = 0
    if (len(Conn) == 0):
        for n in range(0, nnodes + 1):
            c = []
            Conn.append(c)

    # generate mapping used later for resilience calculation
    for i in range(0, (llinks + 1)):
        if (len(Conn) == (nnodes + 1)):
            nodes = d.ENgetlinknodes(i + 1)
            if (nodes[0] == 0):
                Conn[nodes[1]].append(i + 1)
                Conn[nodes[2]].append(i + 1)
                h_degree = max(len(Conn[nodes[0]]), len(Conn[nodes[1]]), h_degree)


    for idx in range(0, nnodes + 1):
        NoConn.append(len(Conn[idx]))
        while (len(Conn[idx]) < h_degree):
            Conn[idx].append(0)
    del Conn[0]
    del NoConn[0]

    return Conn,NoConn
