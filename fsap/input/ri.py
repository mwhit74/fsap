import numpy

def read_input(path):
    
    ifile = open(path, 'r')
    
    read_joint(data)
    read_support(data)
    read_matlprop(data)
    read_sectprop(data)
    read_element(data)
    read_load(data)


def read_joint(data):
    jt = []
    njt = int(data.readline())
    for i in range(njt):
        text = data.readline()
        x_coord = float(text.split(",")[0])
        y_coord = float(text.split(",")[1])
        jt.append((x_coord,y_coord))

    return jt


def read_support(data):
    sup = []
    ns = int(data.readline())
    for i in range(ns):
        text = data.readline()
        jt = int(text.split(",")[0])
        x_dof = int(text.split(",")[1])
        y_dof = int(text.split(",")[2])
        sup.append((jt, x_dof, y_dof))

    return sup


def read_matlprop(data):
    matl = []
    nm = int(data.readline())
    for i in range(nm):
        text = data.readline()
        E = float(text)
        matl.append(E)

    return matl


def read_sectprop(data):
    sp = []
    nsp = int(data.readline())
    for i in range(nsp):
        text = data.readline()
        A = float(text)
        sp.append(A)


def read_element(data):


def read_load(data):
