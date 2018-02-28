import numpy

def read_input(path):
    
    ifile = open(path, 'r')
    
    jt = read_joint(ifile)
    sup = read_support(ifile)
    matl = read_matlprop(ifile)
    sect = read_sectprop(ifile)
    elem = read_element(ifile)
    load = read_load(ifile)

    return jt, sup, matl, sect, elem, load


def read_joint(data):
    jt = []
    njt = int(data.readline())
    for i in range(njt):
        text = data.readline()
        x_coord = float(text.split(",")[0])
        y_coord = float(text.split(",")[1])
        jt.append((x_coord, y_coord)) 
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

    return sp


def read_element(data):
    elem = []
    ne = int(data.readline())
    for i in range(ne):
        text = data.readline()
        beg = int(text.split(",")[0])
        end = int(text.split(",")[1])
        matl = int(text.split(",")[2])
        sp = int(text.split(",")[3])
        elem.append((beg, end, matl,sp))

    return elem


def read_load(data):
    load = []
    nl = int(data.readline())
    for i in range(nl):
        text = data.readline()
        jt = int(text.split(",")[0])
        fx = float(text.split(",")[1])
        fy = float(text.split(",")[2])
        load.append((jt, fx, fy))

    return load
