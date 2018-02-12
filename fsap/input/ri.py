import numpy

def read_input(path):
    
    ifile = open(path, 'r')
    data = ifile.read()
    
    read_joint(data)
    read_support(data)
    read_matlprop(data)
    read_sectprop(data)
    read_element(data)
    read_load(data)

def read_joint(data):
    jt = []
    njt = data.readline()
    for i in range(njt)
        
