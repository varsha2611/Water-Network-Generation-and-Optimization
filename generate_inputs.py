import getopt, sys
import os


def initialize():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:p:g:v", ["input=", "params=", "generator="])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
        sys.exit(2)
    verbose = False
    output_path = None
    params = None
    input_path = None
    job_id = None
    params = {}
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-i", "--input"):
            input_path = a
        elif o in ("-g", "--generator"):
            generator = a
        elif o in ("-p", "--params"):
            try:
                params.update(eval(a.strip()))
            except Exception as inst:
                print('Error parsing parameters!  Given:')
                print(a)
                raise
        else:
            assert False, "unhandled option"
            print(o,a)

    ret = {'params': params, 'input_path': input_path , 'generator':generator}
    return ret


if __name__ == "__main__":
    init_options = initialize()
    input_path = init_options['input_path']
    generator = init_options['generator']
    params = init_options['params']
    output = input_path.replace(".edges", "generated"+generator)
    output = output.replace("input","rescale")
    if input_path == None:
        print("No input network given")
        sys.exit(2)
    file = open('input.txt', "w+")
    for j in range(1,101):
        file.write("python musketeer.py -s "+str(j)+" -f '"+input_path+"' -p \""+str(params)+ "\" -k True -o \""+ output+ str(j)+".edges"+"\"\n")
    file.close()

   




