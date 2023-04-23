
import re
import numpy as np
import string

def read_input(ifn):

    njlc = 0
    ntlc = 0
    ndlc = 0
    ntrlc = 0

    jlcm = np.empty((0,4))
    tlcm = np.empty((0,3))
    dlcm = np.empty((0,4))
    trlcm = np.empty((0,4))
    blcm = np.empty((0,4))

    wob = False
    ofile = ""

    with open(ifn,'r') as ifile:
        for line in ifile:
            line = line.strip()

            if line != "\n":

                if line == "ELEMENT TYPE":
                    elem_type = get_next_line(ifile).lower()

                if line == "CONSTITUTIVE":
                    line = get_next_line(ifile)
                    const = int(get_next_line(ifile))

                if line == "NUM NODES":
                    nn = int(get_next_line(ifile))

                if line == "NUM ELEMS":
                    nelx = int(get_next_line(ifile))

                if line == "NUM MATLS":
                    nm = int(get_next_line(ifile))

                if line == "NUM SECTS":
                    nx = int(get_next_line(ifile))

                if line == "NUM SUPS":
                    ns = int(get_next_line(ifile))

                if line == "NUM LOAD CASES":
                    nlc = int(get_next_line(ifile))

                if line == "NODE COORDS":
                    line = get_next_line(ifile)
                    line = get_next_line(ifile)
                    cols = len(line.split(" "))
                    nc = np.empty([nn,cols])
                    for i in range(nn):
                        nc[i:] = [float(val) for val in line.split(" ")]
                        line = get_next_line(ifile)
                
                if line == "ELEM CNXN":
                    line = get_next_line(ifile)
                    line = get_next_line(ifile)
                    cols = len(line.split(" "))
                    ecm = np.empty((nelx,cols),dtype='int')
                    for i in range(nelx):
                        ecm[i:] = [int(val) for val in line.split(" ")]
                        line = get_next_line(ifile)

                if line == "MATL TYPE":
                    line = get_next_line(ifile)
                    line = get_next_line(ifile)
                    cols = len(line.split(" "))
                    mp = np.empty((nm,cols))
                    for i in range(nm):
                        mp[i:] = [float(val) for val in line.split(" ")]
                        line = get_next_line(ifile)

                if line == "SECTION TYPE":
                    line = get_next_line(ifile)
                    line = get_next_line(ifile)
                    cols = len(line.split(" "))
                    sp = np.empty((nx,cols))
                    for i in range(nx):
                        sp[i:] = [float(val) for val in line.split(" ")]
                        line = get_next_line(ifile)

                if line == "ELEM ASSIGN":
                    line = get_next_line(ifile)
                    line = get_next_line(ifile)
                    cols = len(line.split(" "))
                    emp = np.empty([nelx,cols],dtype='int')
                    for i in range(nelx):
                        emp[i:] = [int(val) for val in line.split(" ")]
                        line = get_next_line(ifile)

                if line == "ELEM PARAM":
                    line = get_next_line(ifile)
                    line = get_next_line(ifile)
                    cols = len(line.split(" "))
                    epa = np.empty((nelx,cols),dtype='int')
                    for i in range(nelx):
                        epa[i:] = [float(val) for val in line.split(" ")]
                        line = get_next_line(ifile)

                if line == "SUPPORTS":
                    line = get_next_line(ifile)
                    line = get_next_line(ifile)
                    cols = len(line.split(" "))
                    sm = np.empty([ns,cols],dtype='int')
                    for i in range(ns):
                        sm[i:] = [int(val) for val in line.split(" ")]
                        line = get_next_line(ifile)

                if 'LOADCASE' in line:
                    lc_name = line.split()[1]
                    line = get_next_line(ifile)
                    while ('LOADCASE' not in line 
                          and line in ["JOINT LOAD","TEMP LOAD","SUP DISP"]
                          and line != ""):
                        load_type = line
                        line = get_next_line(ifile)
                        line = get_next_line(ifile)
                        nl = int(line)
                        line = get_next_line(ifile)

                        if load_type == "JOINT LOAD":
                            njlc = njlc + 1
                            for i in range(nl):
                                line = get_next_line(ifile)
                                jlmt = np.array([np.append([njlc], [float(val) for val in line.split(" ")])])
                                jlcm = np.append(jlcm, jlmt, axis=0)
                        elif load_type == "TEMP LOAD":
                            ntlc = ntlc + 1
                            for i in range(nl):
                                line = get_next_line(ifile)
                                tlmt = np.array([np.append([ntlc], [float(val) for val in line.split(" ")])])
                                tlcm = np.append(tlcm, tlmt, axis=0)
                        elif load_type == "SUP DISP":
                            ndlc = ndlc + 1
                            for i in range(nl):
                                line = get_next_line(ifile)
                                dlmt = np.array([np.append([ndlc], [float(val) for val in line.split(" ")])])
                                dlcm = np.append(dlcm, dlmt, axis=0)

                        line = get_next_line(ifile)

                if "SAVE OUTPUT" == line:
                    line = get_next_line(ifile)
                    line = get_next_line(ifile)
                    wob = line

                if "OUTPATH" == line:
                    line = get_next_line(ifile)
                    line = get_next_line(ifile)
                    ofile = line


    return (elem_type, const, nn, nelx, nm, nx, ns, nc, ecm, mp, sp, 
            emp, epa, sm, njlc, jlcm, ntlc, tlcm, ndlc, dlcm, ntrlc, trlcm, blcm, wob, ofile)

def get_next_line(ifile):

    line = ifile.readline()
   
    #read blank lines 
    while line == "\n":
        line = ifile.readline()

    return line.strip()
