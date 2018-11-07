from epynet import *
import epanettools
from epanettools.epanettools import EPANetSimulation,Node, Link, Pattern
import numpy as np
import time
import random
import warnings
nnodes = 0
llinks = 0
timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())
Pump_energy = 0
#data = [[1 ,76.2 ,0.74],[2, 101.6,0.58],[ 3,152.4,0.41],[4,203.2,0.25],[5, 254,0.15],[6,304.8,0.1],[7,355.6,0.08],[8,406.4,0.06],[9,457.2,0.05],\
 #       [10,508,0.04],[11,609.6,0.03],[12,762,0.025],[13, 914.4, 0.02],[14,1066.8, 0.015],[15,1219.2,0.01],[16,1371.6,0.009],[17,1524,0.008],[18,1625.6,0.007]]

data =[6,8,10,12,14,16,18,20,24,30]
C = 0
nb_run = 0
def SetVariables(d):
    global nnodes
    global llinks
    ret, nnodes = d.ENgetcount(d.EN_NODECOUNT)
    ret, llinks = d.ENgetcount(d.EN_LINKCOUNT)


def Res(pipes, pumps,tanks_diam,tanks_max,tanks_min, d ,hStar,Curves,Conn,NoConn,max_elevation):
    #print ('function start')
    global nnodes
    global llinks
    global C
    global Pump_energy
    tstep = 3600
    global nb_run
    nb_run+=1
    #print (nb_run)
    h = [0] * nnodes
    D = [0] * nnodes
    R = [0] * nnodes
    H = [0] * nnodes
    E = [0] * nnodes
    L = [0] * llinks
    F = [0] * llinks
    Diam = [0] * llinks
    ret, Def_Patterns = d.ENgetcount(d.EN_PATCOUNT)
    idx1 = 0
    idx2 = 0

    #print ('setting curve and diameter')
    # set values for pipe diameter and pump energy
    demand_pattern = []
    #demand_pattern = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]


    for i in xrange(llinks):
        ret, type = d.ENgetlinktype(i + 1)

        if type == 1:
            idx = int(round(pipes[idx1]))
            ret = d.ENsetlinkvalue(i + 1, d.EN_DIAMETER,data[idx])
            idx1+=1

        if type == 2:
            #print("curves",i+1,idx)
            idx = int(round(pumps[idx2]))
            ret= d.ENsetlinkvalue(i + 1, d.EN_STATUS,1)
            #print (x)
            ret = d.ENsetlinkvalue(i + 1, d.EN_CURVE,idx)
            ret, index = d.ENgetlinkvalue(i + 1, d.EN_CURVE)
            #print (index)
            idx2+=1


    i=0
    p = 0

    #print (Def_Patterns)
    for pattern in xrange(Def_Patterns):
        if pattern==0:
                ret, nb = d.ENgetpatternlen(pattern+1)
                for j in xrange(nb):
                    ret, val = d.ENgetpatternvalue(pattern+1, j+1)
                    demand_pattern.append(val)

    idx = 0
    # set value for tanks
    #print ('setting tanks')
    for i in xrange(nnodes):
        ret, type = d.ENgetnodetype(i + 1)
        if type == 2: #Tanks
            ret = d.ENsetnodevalue(i + 1, d.EN_MINLEVEL,0)

    for i in xrange(nnodes):
        ret, type = d.ENgetnodetype(i + 1)
        if type == 2: #Tanks
            elevation = max_elevation+100
            hStar[i] = elevation + tanks_min[idx]
            diameter = tanks_diam[idx]
            min_lvl = tanks_min[idx]
            max_lvl = tanks_max[idx]
            ret = d.ENsetnodevalue(i + 1, d.EN_ELEVATION, elevation)
            ret = d.ENsetnodevalue(i + 1, d.EN_TANKDIAM, diameter)
            if(min_lvl>=0 and min_lvl<max_lvl):
                ret = d.ENsetnodevalue(i + 1, d.EN_TANKLEVEL, min_lvl)
                ret = d.ENsetnodevalue(i + 1, d.EN_MAXLEVEL, max_lvl)
                ret, level = d.ENgetnodevalue(i + 1, d.EN_TANKLEVEL)
                ret = d.ENsetnodevalue(i + 1, d.EN_MINLEVEL, level)
            idx += 1

    # start simulation
    Energy = [0] * llinks
    d.ENopenH()
    d.ENinitH(0)
    idx = 0
    hr= 0
    sim_time = 0
    weighted_sum_res = 0
    sum_weight = 0
    Delta = []
    Penalty = []
    min_h = '+inf'
    valid = True
    #d.ENsaveinpfile('non_feasible'+timeNow()+'.inp')
    while tstep > 0:
        ret, t = d.ENrunH()
        if(ret!=0):
            valid = False
        total_energy = 0
        Numerator = 0
        Denomin_energy = 0
        #print ('values for nodes')
        for i in xrange(nnodes):
            ret,p=d.ENgetnodevalue(i+1, d.EN_PRESSURE )
            h[i] = p
            ret, p = d.ENgetnodevalue(i + 1, d.EN_HEAD)
            H[i] = p
            ret, p = d.ENgetnodevalue(i + 1, d.EN_ELEVATION)
            E[i] = p
            ret, type = d.ENgetnodetype(i+1)

            if type==0: #Junctions
                ret, p = d.ENgetnodevalue(i + 1, d.EN_DEMAND)
                D[i] = p*0.227
                Numerator+=D[i]*(H[i]-hStar[i])*0.3048
                Denomin_energy += D[i] * (H[i]) * 0.3048
                Delta.append(H[i] - hStar[i])
            #print ('values for reservoir')
            if type == 1: #Reservoir
                reser = i+1
                flow = Conn[reser-1][0]
                ret, f_type = d.ENgetlinktype(flow)
                if f_type == 2:
                    linked_nodes = d.ENgetlinknodes(flow)
                    if (linked_nodes[0] == 0):
                        if linked_nodes[1] == reser:
                            next_link = linked_nodes[2]
                        else:
                            next_link = linked_nodes[1]
                        if next_link>0:
                            conn_links = Conn[next_link-1]
                            for link in conn_links:
                                ret, c_type = d.ENgetlinktype(link)
                                if c_type == 1:
                                    flow = link
            #print ('values for tank')
            if type ==2: #Tanks
                #print ('values for tank')
                ret, p =d.ENgetnodevalue(i + 1, d.EN_DEMAND)
                ret, min_level = d.ENgetnodevalue(i + 1, d.EN_MINLEVEL)
                ret, max_level = d.ENgetnodevalue(i + 1, d.EN_MAXLEVEL)
                ret, elevation = d.ENgetnodevalue(i + 1, d.EN_ELEVATION)
                min_head = min_level + elevation
                max_head = max_level + elevation
                if (min_head > H[i]):
                    Delta.append(H[i] - min_head)
                if (max_head < H[i]):
                    Delta.append(max_head - H[i])
                D[i] = p*0.227
                if D[i]>=0:
                    Numerator += D[i] * ((H[i] - hStar[i]) * 0.3048)
                else:
                    Denomin_energy+=(abs(D[i])*H[i]*0.3048)
        #print ('values for links')
        for i in xrange(llinks):
            ret, type = d.ENgetlinktype(i+1)
            ret,p = d.ENgetlinkvalue(i+1, d.EN_FLOW)
            F[i] = p*0.227
            if type == 1:
                ret,p = d.ENgetlinkvalue(i+1, d.EN_DIAMETER)
                Diam[i] = p
                ret,p = d.ENgetlinkvalue(i+1, d.EN_LENGTH)
                L[i] = p
            if type == 2:
                ret,Energy[i] = d.ENgetlinkvalue(i+1, d.EN_ENERGY)
                total_energy += Energy[i]
                Pump_energy = total_energy

        #print ('values for Delta')
        for i in range(0, len(Delta)):
                if Delta[i] <= 0:
                    Penalty.append(Delta[i])
                else:
                    Penalty.append(0)
        unit_wt = 9.81
        energy = (total_energy * 3600)/unit_wt
        Reservoir_Energy = abs(F[flow-1]) * H[reser-1] * 0.3048
        Denomin = np.subtract(Reservoir_Energy, Denomin_energy) + energy

        res = Numerator / Denomin
        weighted_sum_res += res*demand_pattern[hr]
        sum_weight +=demand_pattern[hr]

                #new demand pattern after 3600 seconds
        if(sim_time%3600==0 and sim_time!=0): #startpoint and first hour have same demand pattern
            hr+=1
        ret, tstep = d.ENnextH()
        #print (tstep)
        sim_time += tstep
        idx+=1

    d.ENcloseH()
    #d.clean()
    #print ('simulation end')
    C=0
    #print (min_h)
    if(valid):
        C = abs(np.sum(Penalty))
    else:
        C = abs(np.sum(Penalty))
        if C==0:
            C=100000
        #d.ENsaveinpfile('feasible_solution' + timeNow() + '.inp')


    #print ('setting sum_weight')
    if(sum_weight!=0):
        weighted_res = weighted_sum_res/sum_weight
        #print ('setting sum_weight done')
        #print ('simulation end')
        return weighted_res

    return 0

def Constraint():
    global C
    print("iteration complete")
    return C


def ConstraintTankLevel(d,tanks_max,tanks_min):
    global nnodes
    idx = 0
    for i in xrange(nnodes):
        ret, type = d.ENgetnodetype(i + 1)
        if type == 2: #Tanks
            ret, level = d.ENgetnodevalue(i + 1, d.EN_TANKLEVEL)
            if level < int(round(tanks_min[idx])) or level > int(round(tanks_max[idx])):
                #print("tank level Constraint")
                return 1000
            idx += 1

    return -1000

def Cost(d):

    import numpy as np
    import math

    sum = 0
    pipe_sizes = [6, 8, 10, 12, 14, 16, 18, 20, 24, 30]
    pipe_cost = [137.76,191.552,242.064,314.224,389.664,469.04,554.32,646.816,828.528,1135.208]

    global Pump_energy
    Cpv = 0.10367 #the    present   value    factor
    Cp = 0.12 #the    cost    of    electricity $0.12 / kWh

    NoOfHours = 24

    Nop = 365 *NoOfHours

    OC = Cpv * Nop * Cp * Pump_energy

    sum+=OC

    Tv = [50000 ,100000 ,250000, 500000 ,1000000]  # Tank    Volume
    c = [115000 ,145000 ,325000 ,425000 ,600000]  # Tank    Cost
    for i in xrange(llinks):
        ret, type = d.ENgetlinktype(i + 1)
        if type == 1:
            ret, diam = d.ENgetlinkvalue(i + 1, d.EN_DIAMETER)
            ret, length = d.ENgetlinkvalue(i + 1, d.EN_LENGTH)
            sum += (np.interp(diam, pipe_sizes, pipe_cost) * length)

    #Add cost for tanks
    for i in xrange(nnodes):
        ret, type = d.ENgetnodetype(i + 1)
        if type == 2:
            ret, diameter = d.ENgetnodevalue(i + 1, d.EN_TANKDIAM)
            ret, min_lvl = d.ENgetnodevalue(i + 1, d.EN_MINVOLUME)
            ret, max_lvl = d.ENgetnodevalue(i + 1, d.EN_MAXLEVEL)
            xq = np.pi * 7.48052 * max_lvl * math.pow(diameter,2)
            sum+= np.interp(xq, Tv, c)

    Pump_energy = 0

    return sum

def TankConstraint(tanks_max,tanks_min):
    for idx in range(0,len(tanks_max)):
        if int(round(tanks_min[idx]))-int(round(tanks_max[idx])) >=0:
            return (1000*tanks_min[idx]-tanks_max[idx])

    return -1000

def WriteFeasibleSolution(pipes,pumps,tanks_diam,tanks_max,tanks_min, d ,max_elevation):
    #print ('function start')
    global nnodes
    global llinks
    ret, Def_Patterns = d.ENgetcount(d.EN_PATCOUNT)
    idx1 = 0
    idx2 = 0

    demand_pattern = []



    for i in xrange(llinks):
        ret, type = d.ENgetlinktype(i + 1)

        if type == 1:
            idx = int(round(pipes[idx1]))
            ret = d.ENsetlinkvalue(i + 1, d.EN_DIAMETER,data[idx])
            idx1+=1

        if type == 2:
            #print("curves",i+1,idx)
            idx = int(round(pumps[idx2]))
            ret= d.ENsetlinkvalue(i + 1, d.EN_STATUS,1)
            #print (x)
            ret = d.ENsetlinkvalue(i + 1, d.EN_CURVE,idx)
            ret, index = d.ENgetlinkvalue(i + 1, d.EN_CURVE)
            #print (index)
            idx2+=1


    i=0
    p = 0

    #print (Def_Patterns)
    for pattern in xrange(Def_Patterns):
        if pattern==0:
                ret, nb = d.ENgetpatternlen(pattern+1)
                for j in xrange(nb):
                    ret, val = d.ENgetpatternvalue(pattern+1, j+1)
                    demand_pattern.append(val)

    idx = 0
    # set value for tanks
    for i in xrange(nnodes):
        ret, type = d.ENgetnodetype(i + 1)
        if type == 2: #Tanks
            ret = d.ENsetnodevalue(i + 1, d.EN_MINLEVEL,0)

    for i in xrange(nnodes):
        ret, type = d.ENgetnodetype(i + 1)
        if type == 2: #Tanks
            elevation = max_elevation+100
            diameter = tanks_diam[idx]
            min_lvl = tanks_min[idx]
            max_lvl = tanks_max[idx]
            ret = d.ENsetnodevalue(i + 1, d.EN_ELEVATION, elevation)
            ret = d.ENsetnodevalue(i + 1, d.EN_TANKDIAM, diameter)
            if(min_lvl>=0 and min_lvl<max_lvl):
                ret = d.ENsetnodevalue(i + 1, d.EN_TANKLEVEL, min_lvl)
                ret = d.ENsetnodevalue(i + 1, d.EN_MAXLEVEL, max_lvl)
                ret, level = d.ENgetnodevalue(i + 1, d.EN_TANKLEVEL)
                ret = d.ENsetnodevalue(i + 1, d.EN_MINLEVEL, level)
            idx += 1

    d.ENsaveinpfile('feasible_solution' + timeNow() + '.inp')
