
def read_input(ifn):

    ifile = open(ifn,'r')

    while ifile != EOF:
        line = strip(readline(ifile))

        if line == "STRUCTURE TYPE":
            elem_type = strip(readline(ifile))

        if line == "NUM NODES":
            nn = strip(readline(ifile))

        if line == "NUM ELEMS":
            nelx = strip(readline(ifile))

        if line == "NUM MATLS":
            nm = strip(readline(ifile))

        if line == "NUM SECTS":
            nx = strip(readline(ifile))

        if line == "NUM SUPS":
            ns = strip(readline(ifile))

        if line == "NODE COORDS":
            line = readline(ifile)
            for i in linspace(1,nn):
                line = strip(readline(ifile))
                nc.append(split(line," "))
        
        if line == "ELEM CNXN":
            line = readline(ifile)
            for i in linspace(1,nelx):
                line = strip(readline(ifile))
                ecm.append(split(line," "))

        if line == "SECTION TYPE":
            line = readline(ifile)
            for i in linspace(1,nx):
                line = strip(readline(ifile))
                sp.append(split(line," "))

        if line == "ELEM ASSIGN":
            line = readline(ifile)
            for i in linspace(1,nelx):
                line = strip(readline(ifile))
                emp.append(split(line," "))

        if line == "SUPPORTS":
            line = readline(ifile)
            for i in linspace(1,ns):
                line = strip(readline(ifile))
                sm.append(split(line," "))

        if re.match(line, 'LOADCASE \w+'):
            load_type = strip(readline(ifile))
            line = readline(ifile)
            nl = strip(readline(ifile))
            line = readline(ifile)

            if load_type == "JOINT LOAD":
                njlc = njlc + 1
                for i in linspace(1,nl):
                    line = strip(readline(ifile))
                    jlmt.append(split(line," "))
                jlcm.append(jlmt)
            else if load_type == "TEMP LOAD":
                ntlc = ntlc + 1
                for i in linspace(1,nl):
                    line = strip(readline(ifile))
                    tlmt.append(split(ifile," "))
                tlcm.append(tlmt)

    return (elem_type, nn, nelx, nm, nx, ns, nc, ecm, ap, 
            emp, sm, njlc, jlcm, ntlc, tlcm)
