import getopt, sys
import WDS_services
import algorithms
import  time
import os
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
import multiprocessing

from functools import partial
from contextlib import contextmanager

@contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()

timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())

def initialize():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:p:o:v", ["input=", "params=","output="])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
        sys.exit(2)
    verbose = False
    output_path = None
    params = None
    input_path = None
    params = {}
    print(len(opts))
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-i", "--input"):
            input_path = a
        elif o in ("-p", "--params"):
            try:
                params.update(eval(a.strip()))
            except Exception as inst:
                print('Error parsing parameters!  Given:')
                print(a)
                raise
        elif o in ("-o", "--output"):
            output_path = a
        else:
            assert False, "unhandled option"
            print(o,a)

    ret = {'params': params, 'input_path': input_path, 'output_path': output_path}
    return ret

def generateAndValidate(id, G, generated, network_data,output_path):
    new_network_data = WDS_services.generate_network_data(generated, G, network_data,id)
    WDS_services.write_inp_file(new_network_data, output_path,id)
    return WDS_services.has_solution(output_path,id)

if __name__ == "__main__":
    init_options = initialize()
    input_path = init_options['input_path']
    params = init_options['params']
    output_path = init_options['output_path']
    if input_path == None:
        print("No input network given")
        sys.exit(2)

    G, network_data = WDS_services.read_inp_file(input_path)
    #WDS_services.plot_graph(G, network_data)
    new_G = algorithms.generate_graph(G, params=params, planar=True)
    #print (new_G)

    if output_path == None:
        t_str = timeNow()
        if not os.path.exists('output'):
            os.mkdir('output')
        if not os.path.isdir('output'):
            raise ValueError('Cannot write to directory "output"')
        output_path = 'output/'+os.path.splitext(os.path.basename(input_path))[0] +timeNow()
    iteration = 0
    #new_network_data = WDS_services.generate_network_data(new_G, G, network_data)
    #data_pairs = [[new_G, G, network_data,new_network_data,output_path,str(1)], [new_G, G, network_data,new_network_data,output_path,str(2)], [new_G, G, network_data,new_network_data,output_path,str(3)], \
     #            [new_G, G, network_data,new_network_data,output_path,str(4)],new_G, G, network_data,new_network_data,output_path,str(5)]
    variables = ["","1","2","3","4","5"]
    terminate = False
    while(terminate == False):
        with poolcontext(processes=6) as pool:
            result_list = pool.map(partial(generateAndValidate, G=G, network_data= network_data,output_path=output_path,generated=new_G),variables)
            print (result_list)
        iteration += 1
        if True in result_list:
            terminate = True
            break


