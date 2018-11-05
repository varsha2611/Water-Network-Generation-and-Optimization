import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
import forceatlas2
import random
from random import randrange
from random import uniform
import os
import numpy as np



def write_inp_file(network_data,output_path,id=""):
    if (id != ""):
        output_path = output_path+id+".inp"
    file = open(output_path, "w+")

    #write title
    file.write("[TITLE]\n")
    file.write("Generated synthetic network;\n")

    #write junctions
    file.write("\n[JUNCTIONS]\n")
    file.write(";ID\tElev\tDemand\tPattern\n")
    for junction in network_data['JUNCTIONS']:
        file.write('JU'+str(junction)+"\t")
        file.write(str(network_data['ELEVATIONS'][int(junction)])+"\t")
        if network_data['DEMANDS'][int(junction)] < 0:
            network_data['DEMANDS'][int(junction)] = 0
        file.write(str(network_data['DEMANDS'][int(junction)])+"\t")
        #file.write(str(1) + "\t")
        file.write(";\n")

    # write reservoirs
    file.write("\n[RESERVOIRS]\n")
    file.write(";ID\tHead\tPattern\n")
    for reservoir in network_data['RESERVOIRS']:
        file.write('R'+str(reservoir)+"\t")
        file.write(str(194.00))
        file.write(";\n")

    #write tanks
    file.write("\n[TANKS]\n")
    file.write(";ID\tElevation\tInitLevel\tMinLevel\tMaxLevel\tDiameter\tMinVol\tVolCurve\n")
    for tank in network_data['TANKS']:
        file.write('T'+str(tank)+"\t")
        file.write(str(network_data['ELEVATIONS'][int(tank)])+"\t3.25\t1.00\t10.00\t30.00\t0.00")
        file.write(";\n")

    #write pipes
    file.write("\n[PIPES]\n")
    file.write(";ID\tNode1\tNode2\tLength\tDiameter\tRoughness\tMinorLoss\tStatus\n")
    idx = 1
    #diameters = [3,4,6,8,10,12,14,16,18,20,24,30,36,42,48,54,60,64]
    diameters = [12]
    roughness = [120,130,140]
    #length = 100,3500

    for pipe in network_data['PIPES']:
            file.write("PI" + str(idx)+"\t")
            if pipe[0] in network_data['JUNCTIONS']:
                source = "JU"+str(pipe[0])
            elif pipe[0] in network_data['TANKS']:
                source = "T"+str(pipe[0])
            elif pipe[0] in network_data['RESERVOIRS']:
                source = "R" + str(pipe[0])
            if pipe[1] in network_data['JUNCTIONS']:
                target = "JU"+str(pipe[1])
            elif pipe[1] in network_data['TANKS']:
                target = "T"+str(pipe[1])
            elif pipe[1] in network_data['RESERVOIRS']:
                target = "R" + str(pipe[1])
            file.write(source + "\t" + target + "\t")
            file.write(str(uniform(100,3500))+"\t")
            file.write(str(random.choice(diameters))+"\t")
            file.write(str(random.choice(roughness)) + "\t;\n")
            idx+=1

    # write pumps
    file.write("\n[PUMPS]\n")
    file.write(";ID\tNode1\tNode2\tParameters\n")
    idx=1
    for pump in network_data['PUMPS']:
        file.write("PU" + str(idx) + "\t")
        if pump[0] in network_data['JUNCTIONS']:
            source = "JU"+str(pump[0])
        elif pump[0] in network_data['TANKS']:
            source = "T"+str(pump[0])
        elif pump[0] in network_data['RESERVOIRS']:
            source = "R" + str(pump[0])
        if pump[1] in network_data['JUNCTIONS']:
            target = "JU"+str(pump[1])
        elif pump[1] in network_data['TANKS']:
            target = "T"+str(pump[1])
        elif pump[1] in network_data['RESERVOIRS']:
            target = "R" + str(pump[1])
        file.write(source + "\t" + target + "\t")
        file.write("HEAD 1 \t")
        file.write("PATTERN " + str(idx+1) +"\t;\n")
        idx+=1

    #write coordinates
    file.write("\n[TIMES]\n")
    file.write("Duration\t\t24:00\n")
    file.write("Hydraulic Timestep\t1:00\n")
    file.write("Quality Timestep\t0:05\n")
    file.write("Pattern Timestep\t1:00\n")
    file.write("Pattern Start\t0:00\n")
    file.write("Report Timestep\t1:00\n")
    file.write("Report Start\t0:00\n")
    file.write("Start Clocktime\t12 am\n")
    file.write("Statistic\t\tNONE\n")

    file.write("\n[COORDINATES]\n")
    file.write(";Node\tX - Coord\tY - Coord\n")
    for junction, position in zip(network_data['JUNCTIONS'], network_data['COORDINATES']):
        file.write("JU"+str(junction))
        file.write("\t" + str(network_data['COORDINATES'][junction][0]) + "\t" + str(network_data['COORDINATES'][junction][1]) + "\n")
    for reservoir, position in zip(network_data['RESERVOIRS'], network_data['COORDINATES']):
        file.write("R"+str(reservoir))
        file.write("\t" + str(network_data['COORDINATES'][reservoir][0]) + "\t" + str(network_data['COORDINATES'][reservoir][1]) + "\n")
    for tank, position in zip(network_data['TANKS'], network_data['COORDINATES']):
        file.write("T"+str(tank))
        file.write("\t" + str(network_data['COORDINATES'][tank][0]) + "\t" + str(network_data['COORDINATES'][tank][1]) + "\n")

    file.write("\n[VALVES]\n")
    file.write(";ID\tNode1\tNode2\tDiameter\tType\tSetting\tMinorLoss\n")
    file.write("\n[PATTERNS]\n")
    file.write(";ID\tMultipliers\n")
    pattern = 1
    file.write(str(pattern) + "\t1.0000\t1.0000\t1.0000\t0.9000\t0.9000\t0.9000\n")
    file.write(str(pattern) + "\t.7000\t.7000\t0.7000\t0.6000\t0.6000\t0.6000\n")
    file.write(str(pattern) + "\t1.2000\t1.2000\t1.2000\t1.3000\t1.3000\t1.3000\n")
    file.write(str(pattern) + "\t1.2000\t1.2000\t1.2000\t1.1000\t1.1000\t1.1000\n")
    pattern+= 1
    for pump in network_data['PUMPS']:
        for hr in xrange(0,4):
            file.write(str(pattern)+"\t1.0000\t1.0000\t1.0000\t1.0000\t1.0000\t1.0000\n")
        pattern+=1
    file.write("\n[TAGS]\n")
    file.write("\n[DEMANDS]\n")
    file.write("\n[STATUS]\n")
    file.write("\n[CURVES]\n")
    file.write(";ID \tX-Value\t Y-Value\n")
    file.write("1\t10\t280\n")
    file.write("\n[CONTROLS]\n")
    file.write("\n[RULES]\n")
    file.write("\n[ENERGY]\n")

    file.close()
    return output_path

def read_inp_file(filename):
    G = nx.Graph()
    file = open(filename, 'r')
    network_data = {}
    elevation = []
    reservoirs = []
    tanks = []
    pumps = []
    junctions = []
    pipes = []
    roughness = []
    lengths = []
    demands = []
    for fulline in file:
        # Network Components
        key = fulline.split()
        if(len(key)>0):
            line = key[0]
        if line == "[TITLE]":
           current_tag = "TITLE"
        elif line == "[PIPES]":
            current_tag = "PIPES"
        elif line == "[JUNCTIONS]":
            current_tag = "JUNCTIONS"
        elif line == "[RESERVOIRS]":
            current_tag = "RESERVOIRS"
        elif line == "[TANKS]":
            current_tag = "TANKS"
        elif line == "[PIPES]":
            current_tag = "PIPES"
        elif line == "[PUMPS]":
            current_tag = "PUMPS"
        elif line == "[VALVES]":
            current_tag = "VALVES"
        elif line == "[EMITTERS]":
            current_tag = "EMITTERS"
        elif line == "[CURVES]":
            current_tag = "CURVES"
        elif line == "[PATTERNS]":
            current_tag = "PATTERNS"
        elif line == "[ENERGY]":
            current_tag = "ENERGY"
        elif line == "[STATUS]":
            current_tag = "STATUS"
        elif line == "[CONTROLS]":
            current_tag = "CONTROLS"
        elif line == "[RULES]":
            current_tag = "RULES"
        elif line == "[DEMANDS]":
            current_tag = "DEMANDS"
        elif line == "[QUALITY]":
            current_tag = "QUALITY"
        elif line == "[REACTIONS]":
            current_tag = "REACTIONS"
        elif line == "[SOURCES]":
            current_tag = "SOURCES"
        elif line == "[MIXING]":
            current_tag = "MIXING"
        elif line == "[OPTIONS]":
            current_tag = "OPTIONS"
        elif line == "[TIMES]":
            current_tag = "TIMES"
        elif line == "[REPORT]":
            current_tag = "REPORT"
        elif line == "[COORDINATES]":
            current_tag = "COORDINATES"
        elif line == "[TAGS]":
            current_tag = "TAGS"
        elif line[0] == ";" or fulline.isspace():
            continue
        elif line[0] == "[":
            current_tag = "INVALID"
        elif current_tag == "PIPES" or current_tag == "VALVES":
            line = fulline.split()
            pipes.append(line[0])
            lengths.append(line[3])
            roughness.append(line[5])
            G.add_edge(str(line[1]),str(line[2]))
        elif current_tag == "PUMPS":
            line = fulline.split()
            pumps.append(line[0])
            G.add_edge(str(line[1]), str(line[2]))
        elif current_tag == "JUNCTIONS":
            line = fulline.split()
            elevation.append(line[1])
            demands.append(line[2])
            junctions.append(str(line[0]))
        elif current_tag == "RESERVOIRS":
            line = fulline.split()
            reservoirs.append(str(line[0]))
        elif current_tag == "TANKS":
            line = fulline.split()
            tanks.append(str(line[0]))
            elevation.append(line[1])
        else:
            current_tag ="NONE"

    network_data['JUNCTIONS'] = junctions
    network_data['PUMPS'] = pumps
    network_data['RESERVOIRS'] = reservoirs
    network_data['TANKS'] = tanks
    network_data['PIPES'] = pipes
    network_data['ELEVATIONS'] = elevation
    network_data['ROUGHNESS'] = roughness
    network_data['LENGTH'] = lengths
    network_data['DEMANDS'] = demands
    print(len(G.edges()),len(G.nodes()))
    file.close()
    return G,network_data

def generate_network(G,output_path,network_data):
    initial_pos =  nx.graphviz_layout(G,prog='neato')
    pos = forceatlas2.forceatlas2_networkx_layout(G, pos = initial_pos, niter=100,gravity=.12,strongGravityMode=True,scalingRatio = 5.0) # Optionally specify iteration count
    probable_reservoirs = []
    min_x = min_y = float('inf')
    max_x = max_y = float('-inf')
    for node in pos:
        if(min_x > pos[node][0]) and (G.degree(node) == 1):
            min_x = pos[node][0]
            probable_reservoirs.append(node)
        if(max_x < pos[node][0]) and (G.degree(node) == 1):
            max_x = pos[node][0]
            probable_reservoirs.append(node)
        if (min_y > pos[node][1]) and (G.degree(node) == 1):
            min_y = pos[node][1]
            probable_reservoirs.append(node)
        if (max_y < pos[node][1]) and (G.degree(node) == 1):
           max_y = pos[node][1]
           probable_reservoirs.append(node)

    random_index = randrange(int(len(probable_reservoirs)/2), len(probable_reservoirs))
    reservoir = probable_reservoirs[random_index]
    node_colors = ['blue' if node == reservoir else 'red' for node in G.nodes()]
    nx.draw(G, pos, with_labels=False, node_color=node_colors, node_size=50)
    write_inp_file(network_data,output_path)
    #G_new , new_data = read_inp_file(output_path)
    #get_data_for_optimization(G_new, new_data)

def assign_reservoirs(G, network_data, new_network_data):

    if network_data.has_key('RESERVOIRS'):
        nb_of_reservoirs = len(network_data['RESERVOIRS'])
        probable_reservoirs = []
        reservoirs = []
        pos = new_network_data['COORDINATES']
        min_x = min_y = float('inf')
        max_x = max_y = float('-inf')
        for node in pos:
            if (min_x > pos[node][0]) and (G.degree(node) == 1):
                min_x = pos[node][0]
                probable_reservoirs.append(node)
            if (max_x < pos[node][0]) and (G.degree(node) == 1):
                max_x = pos[node][0]
                probable_reservoirs.append(node)
            if (min_y > pos[node][1]) and (G.degree(node) == 1):
                min_y = pos[node][1]
                probable_reservoirs.append(node)
            if (max_y < pos[node][1]) and (G.degree(node) == 1):
                max_y = pos[node][1]
                probable_reservoirs.append(node)

        random_index = random.sample(range(int(len(probable_reservoirs) / 2), len(probable_reservoirs)),nb_of_reservoirs)
        b = 1
        for index in random_index:
            probability = random.uniform(0, b)
            if probability >0.5 or len(reservoirs) == 0 :
                reservoirs.append(probable_reservoirs[index])
                b-=0.1
            else:
                b+=0.1
        new_network_data['RESERVOIRS'] = reservoirs

    return new_network_data

def write_metis_format(G,output_path):
    file = open(output_path, "w+")
    file.write(str(len(G.nodes()))+"\t"+ str(len(G.edges()))+"\n")
    for node in G.nodes():
        #file.write(str(node) + "\t")
        neighbors = G.neighbors(node)
        for neighbor in neighbors:
            file.write(str(neighbor)+"\t")
        file.write("\n")
    file.close()

def assign_tanks_and_pumps(new_G, G, network_data, new_network_data,id=""):
    if network_data.has_key('TANKS'):
        rescale = int(len(new_G.nodes())/len(G.nodes()))
        new_network_data['TANKS'] = []
        new_network_data['PUMPS'] = []

        nb_of_tanks = len(network_data['TANKS'])*rescale
        if nb_of_tanks == 0:
            nb_of_tanks = 4

        metis_format_file = 'temp_graph'+id+'.graph'
        write_metis_format(new_G,metis_format_file)
        seed = random.randint(0,len(new_G.nodes()))
        output_file_name = 'partitions'+id
        import subprocess
        subprocess.call(["/home/varsha/Documents/softwares/KaHIP-master/deploy/kaffpaE", metis_format_file, '--k',str(nb_of_tanks),'--seed', str(seed),'--imbalance','10','--preconfiguration','strong','--output_filename',output_file_name])
        Partitions = read_partition(output_file_name)
        os.remove(output_file_name)
        b=1.5
        for partition in Partitions:
            probability = random.uniform(0, b)
            if probability > 0.5:
                for node in Partitions[partition]:
                    if new_G.degree(node) == 1:
                        new_network_data['TANKS'].append(node)
                        b -= 0.1
                        break
                    else:
                        continue

            else:
                b += 0.1


        reservoir = random.choice(new_network_data['RESERVOIRS'])
        b=1.5
        for tank in new_network_data['TANKS']:
            path = nx.shortest_path(new_G, source=reservoir, target=tank, weight=None)
            if reservoir in path and tank in path:
                path.remove(reservoir)
                path.remove(tank)
                if len(path)==0:
                    continue
                index = random.choice(range(0, len(path)))
                probability = random.uniform(0, b)
                source = path[index]
                if index == (len(path)-1):
                    target = tank
                else:
                    target = path[index+1]
                if probability > 0.5:
                    new_network_data['PUMPS'].append((source, target))
                    b -= 0.1
                else:
                    b += 0.1

        return new_network_data


def generate_network_data(new_G,G,network_data,id=""):
    labeled_new_G = nx.convert_node_labels_to_integers(new_G, first_label=1, ordering='default', label_attribute=None)
    new_network_data = {}

    # Assign coordinates
    initial_pos = nx.graphviz_layout(labeled_new_G, prog='neato')
    pos = forceatlas2.forceatlas2_networkx_layout(labeled_new_G, pos=initial_pos, niter=100, gravity=0.12, strongGravityMode=True,scalingRatio=5.0)

    new_network_data['COORDINATES'] = pos

    #Assign reservoirs
    new_network_data = assign_reservoirs(labeled_new_G,network_data,new_network_data)

    # Assign tanks
    new_network_data = assign_tanks_and_pumps(labeled_new_G,G,network_data,new_network_data,id)

    #Assign junctions
    new_network_data = assign_junctions(labeled_new_G,network_data,new_network_data)

    #Assign elevation
    new_network_data = assign_elevation(labeled_new_G,network_data,new_network_data)

    #Assign pipes
    new_network_data = assign_pipes(labeled_new_G,network_data,new_network_data)

    #Assign demand
    new_network_data = assign_demand(labeled_new_G, network_data, new_network_data)

    return new_network_data

def plot_graph(G,network_data):
    initial_pos = nx.graphviz_layout(G, prog='neato')
    pos = forceatlas2.forceatlas2_networkx_layout(G, pos=initial_pos, niter=100, gravity=0.12,strongGravityMode=True, scalingRatio=5.0)
    blue_nodes = []
    red_nodes =[]
    black_nodes =[]
    for node in network_data['RESERVOIRS']:
        print ('reservoir',node)
        blue_nodes.append(str(node))
    for node in network_data['TANKS']:
        print('tank', node)
        red_nodes.append(str(node))
    for node in network_data['JUNCTIONS']:
        print('junction', node)
        black_nodes.append(str(node))

    nx.draw_networkx_nodes(G, pos, nodelist=blue_nodes, node_color='b', node_size=20)
    nx.draw_networkx_nodes(G, pos, nodelist=red_nodes, node_color='r', node_size=20)
    nx.draw_networkx_nodes(G, pos, nodelist=black_nodes, node_color='black', node_size=5)
    #nx.draw_networkx_nodes(G, pos, nodelist=other_nodes, node_color='black', node_size=5)
    nx.draw_networkx_edges(G, pos)
    # nx.draw(G, pos, with_labels=False, node_color=node_colors, node_size=50)
    plt.show()

def read_partition(file_path):
    Partitions = {}
    file = open(file_path)
    node_number = 1
    for line in file:
        values = line.split()
        partition_nb = values[0]
        partition = []
        if Partitions.has_key(partition_nb):
            Partitions[partition_nb].append(node_number)
        else:
            partition.append(node_number)
            Partitions[partition_nb] = partition
        node_number +=1

    return Partitions

def assign_junctions(new_G,network_data,new_network_data):
    junctions = []
    for node in new_G.nodes():
        if node not in new_network_data['TANKS'] and node not in new_network_data['RESERVOIRS']:
            junctions.append(node)

    new_network_data['JUNCTIONS'] = junctions

    return new_network_data


def assign_elevation(new_G, network_data, new_network_data):
    original_distribution = network_data['ELEVATIONS']
    new_distribution = []

    #initialize with zero
    for node in new_G.nodes():
        new_distribution.append(0)

    #neglect zero position add additional indexing
    new_distribution.append(0)

    #Assign initial random values
    for node in new_G.nodes():
        new_distribution[int(node)] = original_distribution[((int(node))%(len(original_distribution)))]

    #smoothing process
    iterations = 3
    smoothened_values = smoothen_values(new_G,new_distribution,iterations)
    #update Elevation for tanks
    for idx in range(1,len(smoothened_values)):
        if idx in new_network_data['TANKS']:
            neighbors = new_G.neighbors(idx)
            sum = 0
            for neighbor in neighbors:
                sum += float(smoothened_values[neighbor])
            smoothened_values[idx] = sum+random.uniform(50, 70)
        idx += 1

    new_network_data['ELEVATIONS'] = smoothened_values
    return new_network_data

def assign_demand(new_G, network_data, new_network_data):
    original_distribution = network_data['DEMANDS']
    new_distribution = []

    # initialize with zero
    for node in new_G.nodes():
        new_distribution.append(0)

    # neglect zero position add additional indexing
    new_distribution.append(0)

    #Assign initial random values
    for node in new_G.nodes():
        new_distribution[int(node)] = original_distribution[int(node)%len(original_distribution)]

    #smoothing process
    iterations = 3
    new_network_data['DEMANDS'] = smoothen_values(new_G,new_distribution,iterations)

    return new_network_data

def assign_pipes(new_G, network_data, new_network_data):
    new_network_data['PIPES'] = []
    for edge in new_G.edges():
        if (edge[0],edge[1]) not in new_network_data['PUMPS'] and (edge[1],edge[0]) not in new_network_data['PUMPS']:
            new_network_data['PIPES'].append((edge[0],edge[1]))

    return new_network_data

def smoothen_values(G,distribution,iterations):
    itr = 0
    while itr < iterations:
        idx = 1
        while idx<len(distribution):
            neighbors = G.neighbors(idx)
            sum = 0
            for neighbor in neighbors:
                sum += float(distribution[neighbor])
            distribution[idx]=float(sum/len(neighbors))
            idx+=1
        itr+=1

    return distribution


def has_solution(input_network,id=""):
    import sys
    from pathlib import Path
    from platypus import NSGAII, Problem, Integer, Real
    import Functions
    import Settings
    class my_mo_problem(Problem):
        et, hStar, o_curves, n_curves, nbOfPipes, nbOfPumps, nbOfTanks, Conn, NoConn, max_elevation = Settings.SetValues(input_network,id=id)
        Functions.SetVariables(et)

        def __init__(self):
            super(my_mo_problem, self).__init__(self.nbOfPipes + 25 * self.nbOfPumps + 3 * self.nbOfTanks, 2, 1)
            self.types[:] = [Real(0, 9)] * self.nbOfPipes + [Real(0, self.n_curves - 1)] * self.nbOfPumps + [Real(25, 100)] * self.nbOfTanks + [Real(25, 40)] * self.nbOfTanks + [Real(9, 10)] * self.nbOfTanks + [Real(0, 1)] * (24 * self.nbOfPumps)
            self.constraints[:] = "<=0"
            self.directions[:] = Problem.MINIMIZE
            # self.function = mixed_type

        def evaluate(self, solution):
            pipes = solution.variables[0:self.nbOfPipes]  # diameter of pipe
            pumps = solution.variables[self.nbOfPipes:self.nbOfPipes + self.nbOfPumps]  # curve of pumps
            tanks_diam = solution.variables[
                         self.nbOfPipes + self.nbOfPumps:self.nbOfPipes + self.nbOfPumps + self.nbOfTanks]  # diameter of tank
            tanks_max = solution.variables[
                        self.nbOfPipes + self.nbOfPumps + self.nbOfTanks:self.nbOfPipes + self.nbOfPumps + 2 * self.nbOfTanks]  # max level of tank
            tanks_min = solution.variables[
                        self.nbOfPipes + self.nbOfPumps + 2 * self.nbOfTanks:self.nbOfPipes + self.nbOfPumps + 3 * self.nbOfTanks]  # min level of tank
            patterns = solution.variables[
                       self.nbOfPipes + self.nbOfPumps + 3 * self.nbOfTanks:self.nbOfPipes + self.nbOfPumps + 3 * self.nbOfTanks + 24 * self.nbOfPumps]
            solution.objectives[:] = [-Functions.Res(pipes, patterns, pumps, tanks_diam, tanks_max, tanks_min, self.et, self.hStar,self.n_curves, self.Conn, self.NoConn, self.max_elevation),Functions.Cost(patterns, self.et)]
            solution.constraints[:] = [Functions.Constraint()]

    algorithm = NSGAII(my_mo_problem())
    algorithm.run(1000)

    feasible_solutions = [s for s in algorithm.result if s.feasible]
    if(len(feasible_solutions)>0):
        nondominated_solutions = nondominated(feasible_solutions)

        sln = 1
        file = open('solution'+id+'.txt', "w+")
        for solution in nondominated_solutions:

            pipes = solution.variables[0:nbOfPipes]  # diameter of pipe
            pumps = solution.variables[nbOfPipes:nbOfPipes + nbOfPumps]  # curve of pumps
            tanks_diam = solution.variables[
                         nbOfPipes + nbOfPumps:nbOfPipes + nbOfPumps + nbOfTanks]  # diameter of tank
            tanks_max = solution.variables[
                        nbOfPipes + nbOfPumps + nbOfTanks:nbOfPipes + nbOfPumps + 2 * nbOfTanks]  # max level of tank
            tanks_min = solution.variables[
                        nbOfPipes + nbOfPumps + 2 * nbOfTanks:nbOfPipes + nbOfPumps + 3 * nbOfTanks]  # min level of tank
            patterns = solution.variables[
                       nbOfPipes + nbOfPumps + 3 * nbOfTanks:nbOfPipes + nbOfPumps + 3 * nbOfTanks + 24 * nbOfPumps]
            Functions.WriteFeasibleSolution(pipes, patterns, pumps, tanks_diam, tanks_max, tanks_min, et, max_elevation)
            file.write("solution index " + str(sln))
            file.write(
                "\nResilience: " + str(-solution.objectives[0]) + "\tCost: " + str(solution.objectives[1]) + "\n")
            file.write("\npipes\n")
            for i in range(0, len(pipes)):
                pipes[i] = int(round(pipes[i]))

            file.write(str(pipes))

            file.write("\npumps\n")
            for i in range(0, len(pumps)):
                pumps[i] = int(round(pumps[i]))

            file.write(str(pumps))

            for i in range(0, len(patterns)):
                patterns[i] = int(round(patterns[i]))

            for i in range(0, len(tanks_max)):
                tanks_max[i] = tanks_max[i]
            file.write("\ntank diameters\n")
            file.write(str(tanks_diam))

            file.write("\ntank maximum level\n")
            file.write(str(tanks_max))

            for i in range(0, len(tanks_min)):
                tanks_min[i] = tanks_min[i]
            file.write("\ntank minimum level\n")
            file.write(str(tanks_min))

            file.write("\npatterns\n")
            p = 0
            idx = 0
            for pattern in patterns:
                if (idx >= 1 and idx <= len(pumps)):
                    file.write(str(patterns[p:p + 24]) + "\n")
                    p += 23
                idx += 1
            sln += 1
            file.write("\n\n")


def generate_PlatypusProblem(input_network):
    problem = open("WDS_Scripts/MyCustomProblem.py", "w+")
    problem.write("from platypus import NSGAII, Problem, Integer,Real\n")
    problem.write("import sys\n")
    problem.write("sys.path.insert(0, '/home/varsha/Documents/MyCode/Water Network/Optimization/')\n")
    problem.write("import Functions\n")
    problem.write("import Settings\n")
    problem.write("class my_mo_problem(Problem):\n")
    problem.write("\tet, hStar, o_curves, n_curves, nbOfPipes,nbOfPumps,nbOfTanks,Conn,NoConn,max_elevation = Settings.SetValues('"+input_network+"')\n")
    problem.write("\tFunctions.SetVariables(et)\n")
    problem.write("\tdef __init__(self):\n")
    problem.write("\t\tsuper(my_mo_problem, self).__init__(self.nbOfPipes + 25*self.nbOfPumps + 3*self.nbOfTanks ,2, 1)\n")
    problem.write("\t\tself.types[:] = [Real(0, 9)]*self.nbOfPipes + [Real(0,self.n_curves-1)]*self.nbOfPumps + [Real(25,100)]*self.nbOfTanks + [Real(25,40)]*self.nbOfTanks+ [Real(9,10)]*self.nbOfTanks+[Real(0,1)]*(24*self.nbOfPumps)\n")
    problem.write("\t\tself.constraints[:] = \"<=0\"\n")
    problem.write("\t\tself.directions[:] = Problem.MINIMIZE\n")
    problem.write("\tdef evaluate(self, solution):\n")
    problem.write("\t\tpipes = solution.variables[0:self.nbOfPipes] #diameter of pipe\n")
    problem.write("\t\tpumps = solution.variables[self.nbOfPipes:self.nbOfPipes+self.nbOfPumps] #curve of pumps\n")
    problem.write("\t\ttanks_diam = solution.variables[self.nbOfPipes+self.nbOfPumps:self.nbOfPipes+self.nbOfPumps+self.nbOfTanks] #diameter of tank\n")
    problem.write("\t\ttanks_max =solution.variables[self.nbOfPipes+self.nbOfPumps+self.nbOfTanks:self.nbOfPipes+self.nbOfPumps+2*self.nbOfTanks] # max level of tank\n")
    problem.write("\t\ttanks_min = solution.variables[self.nbOfPipes+self.nbOfPumps+2*self.nbOfTanks:self.nbOfPipes+self.nbOfPumps+3*self.nbOfTanks] #min level of tank\n")
    problem.write("\t\tpatterns = solution.variables[self.nbOfPipes + self.nbOfPumps + 3 * self.nbOfTanks:self.nbOfPipes + self.nbOfPumps + 3 * self.nbOfTanks + 24 * self.nbOfPumps]\n")
    problem.write("\t\tsolution.objectives[:] = [-Functions.Res(pipes,patterns,pumps,tanks_diam,tanks_max,tanks_min,self.et,self.hStar,self.n_curves,self.Conn,self.NoConn,self.max_elevation),Functions.Cost(patterns,self.et)]\n")
    problem.write("\t\tsolution.constraints[:] = Functions.Constraint()\n")

def generate_PlatypusProblem2(input_network):
    problem = open('/home/varsha/Documents/MyCode/Water Network/Optimization/MyProblem.py', "w+")
    problem.write("from platypus import NSGAII, Problem, Integer\n")
    problem.write("import Constraint\n")
    problem.write("import Cost\n")
    problem.write("import Res\n")
    problem.write("from epanettools.epanettools import EPANetSimulation\n")
    problem.write("class my_mo_problem(Problem):\n")
    problem.write("\td = EPANetSimulation('"+input_network+"')\n")
    problem.write("\tret, nlinks = d.ENgetcount(d.EN_LINKCOUNT)\n")
    problem.write("\tdef __init__(self):\n")
    problem.write("\t\tsuper(my_mo_problem, self).__init__(self.nlinks, 2, 2)\n")
    problem.write("\t\tself.types[:] = [Integer(0, 16)]*self.nlinks\n")
    problem.write("\t\tself.constraints[:] = \"<=0\"\n")
    problem.write("\t\tself.directions[:] = Problem.MINIMIZE\n")
    problem.write("\tdef evaluate(self, solution):\n")
    problem.write("\t\ty = solution.variables\n")
    problem.write("\t\tsolution.objectives[:] = [Res.Res(y,self.d),Cost.Cost(y,self.d)]\n")
    problem.write("\t\tsolution.constraints[:] = Constraint.Constraint(y,self.d)\n")


