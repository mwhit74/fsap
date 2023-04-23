import numpy as np
from numpy import matmul,dot
from fsap.fe_processor import elementdata, calc_dofs
from fsap.fe_processor import (bar221i, bar321i, tri231i, tri263i, quad244i,
                         quad289i, tet344i, cube388i)
import pathlib

def post(inp,K,Kev,F,u,ifile):

    (elem_type, const, nn, nelx, nm, nx, ns, nc, ecm, mp, sp, 
            epm, epa, sm, njlc, jlcm, ntlc, tlcm, ndlc, dlcm, ntrlc, trlcm, blcm, wob, ofile) = inp

    strains, stresses = calc_es(elem_type,const,nelx,nc,ecm,mp,sp,epm,epa,u)

    reacs = calc_reacs(K,F,u)

    write_output_file(ofile, ifile, elem_type, nelx, nn, u, strains, stresses, reacs, wob)

    return strains,stresses,reacs

    
    
def calc_es(elem_type,const,nelx,nc,ecm,mp,sp,epm,epa,u):

    """Calculate the strains and stresses

    Params:
        elem_type (str): unique type
        nelx (int): number of global elements
        nc (ndarray): nodal coordinate matrix
        ecm (ndarray): element connection matrix
        mp (ndarray): material property matrix
        sp (ndarray): section properties matrix
        epm (ndarray): assigned element property matrix
        K (ndarray): global stiffness matrix w/ BCs applied
        F (ndarray): global load vector w/ BCs applied
        u (ndarray): nodal displacment vector

    Returns:
        stresses (ndarray): element stresses
        reacs (ndarray): nodal reactions for restrained nodes
    """

    ssd = elementdata[elem_type][6]
    nip = elementdata[elem_type][2]
    qt_dict = elementdata[elem_type][4]
    qt_id = elementdata[elem_type][5]

    qt = qt_dict[qt_id]
    ips = qt[0]

    stresses = np.empty((nelx,nip,ssd))
    strains = np.empty((nelx,nip,ssd))

    for i in range(nelx):
        nodes = ecm[i,1:]
        elem_fnc = elementdata[elem_type][3]
        dofs = calc_dofs(elem_type, nodes)

        ue = u[dofs]

        for j in range(nip):

            ip = ips[j]

            S,D,B = elem_fnc(i,const,nc,ecm,mp,sp,epm,epa,ip)

            if elem_type in ["bar221i","bar321i"]:

                strains[i,j] = dot(B,ue)
                stresses[i,j] = D*dot(B,ue)

            elif elem_type in ["tri231i","quad244i","quad289i"]:

                strains[i,j] = matmul(B,ue)
                stresses[i,j] = matmul(matmul(D,B),ue)

            elif elem_type in ["tet344i","cube388i"]:

                s = matmul(B,ue)
                st = matmul(matmul(D,B),ue)

                #s11, s12, s13, s22, s23, s33
                strains = np.array([s[0,0],s[0,1],s[0,2],s[1,1],s[1,2],s[2,2]])
                stresses = np.array([st[0,0],st[0,1],st[0,2],st[1,1],st[1,2],st[2,2]])

    return strains, stresses

def calc_reacs(K,F,u):

    reacs = matmul(K,u)-F

    return reacs

def write_output_file(ofile, ifile, elem_type, nelx, nn, 
                        disps, strains, stresses, reacs, wob):

    if wob == True:
        if ifile != "" and ofile == "":
            ipath = pathlib.Path(ifile)
            ofile = pathlib.Path.joinpath(ipath.parent,ipath.stem+"_out.txt")
        if ofile == "":
            raise ValueError("No output file path provided.")

        with open(ofile, 'w') as ofl:
            echo_input_file(ofl, ifile)

        with open(ofile, 'a') as ofl:

            ofl.write("\n\n\n-----------OUTPUT------------\n")

            ofl.write("\n\nNodal Displacements\n")



            if elem_type in ["bar221i","tri231i","quad244i","quad289i"]:
                disp_header = f"{'Node':^10s}{'dx':^10s}{'dy':^10s}\n"
            elif elem_type in ["bar231i","tet344i","cube388i"]:
                disp_header = f"{'Node':^10s}{'dx':^10s}{'dy':^10s}{'dz':^10s}\n"
            elif elem_type in ["beam"]:
                disp_header = (f"{'Node':^10s}{'dx':^10s}{'dy':^10s}{'dz':^10s}"
                                            f"{'rx':^10s}{'ry':^10s}{'rz':^10s}\n")
        
            ofl.write(disp_header)

            j = 0
            for i in range(nn): 
                line = ""
                if elem_type in ["bar221i","tri231i","quad244i","quad289i"]:
                    line = f"{i+1:^10d}{disps[j]:^10.6f}{disps[j+1]:^10.6f}"
                    j = j+2    
                elif elem_type in ["bar321i","tet344i","cube388i"]:
                    line = f"{i+1:^10d}{disps[i]:^10.6f}{disps[i+1]:^10.6f}{disps[i+2]:^10.6f}"
                    j = j+3
                elif elem_type in ["beam"]:
                    line = (f"{i+1:^10d}{disps[j]:^10.6f}{disps[j+1]:^10.6f}{disps[j+2]:^10.6f}"
                                     f"{disps[j+3]:^10.6f}{disps[j+4]:^10.6f}{disps[j+5]:^10.6f}")
                    j = j+6

                ofl.write(line+"\n")



            nip = elementdata[elem_type][2]
            ssd = elementdata[elem_type][6]



            ofl.write("\n\nElement Strains\n")

            if elem_type in ["bar221i","bar321i"]:
                rs_header = f"{'Element':^10s}{'IP':^10s}{'eps_xx':^10s}"
            elif elem_type in ["tri231i","quad244i","quad289i"]:
                rs_header = (f"{'Element':^10s}{'IP':^10s}"
                            f"{'eps_xx':^10s}{'eps_yy':^10s}{'tau_xy':^10s}")
            elif elem_type in ["tet344i","cube388i"]:
                rs_header = (f"{'Element':^10s}{'IP':^10s}"
                            f"{'eps_xx':^10s}{'eps_yy':^10s}{'eps_zz':^10s}"
                            f"{'tau_xy':^10s}{'tau_xz':^10s}{'tau_yz':^10s}")
            elif elem_type in ["beam"]:
                rs_header = f"{'Element':^10s}{'IP':^10s}"

            ofl.write(rs_header+"\n")

            for i in range(nelx):
                line = f"{i+1:^10d}"
                for j in range(nip):
                    t = 'IP ' + str(j)
                    line = line + f"{t:^10s}"
                    for k in range(ssd):
                        line = line + f"{strains[i,j,k]:^10.6f}"
                ofl.write(line + "\n")



            ofl.write("\n\nRaw Element Stresses\n")
            
            if elem_type in ["bar221i","bar321i"]:
                rs_header = f"{'Element':^10s}{'IP':^10s}{'sigma_xx':^10s}"
            elif elem_type in ["tri231i","quad244i","quad289i"]:
                rs_header = (f"{'Element':^10s}{'IP':^10s}"
                            f"{'sigma_xx':^10s}{'sigma_yy':^10s}{'tau_xy':^10s}")
            elif elem_type in ["tet344i","cube388i"]:
                rs_header = (f"{'Element':^10s}{'IP':^10s}"
                            f"{'sigma_xx':^10s}{'sigma_yy':^10s}{'sigma_zz':^10s}"
                            f"{'tau_xy':^10s}{'tau_xz':^10s}{'tau_yz':^10s}")
            elif elem_type in ["beam"]:
                rs_header = f"{'Element':^10s}{'IP':^10s}"

            ofl.write(rs_header+"\n")

            for i in range(nelx):
                line = f"{i+1:^10d}"
                for j in range(nip):
                    t = 'IP ' + str(j)
                    line = line + f"{t:^10s}"
                    for k in range(ssd):
                        line = line + f"{stresses[i,j,k]:^10.2f}"
                ofl.write(line + "\n")
                


            ofl.write("\n\nNodal Reactions\n")
                    
            if elem_type in ["bar221i","tri231i","quad244i","quad289i"]:
                reac_header = f"{'Node':^10s}{'Rx':^10s}{'Ry':^10s}\n"
            elif elem_type in ["tet344i","cube388i"]:
                reac_header = f"{'Node':^10s}{'Rx':^10s}{'Ry':^10s}{'Rz':^10s}\n"
            elif elem_type in ["beam"]:
                reac_header = (f"{'Node':^10s}{'Rx':^10s}{'Ry':^10s}{'Rz':^10s}"
                                            f"{'Mx':^10s}{'My':^10s}{'Mz':^10s}\n")
        
            ofl.write(reac_header)

            j = 0
            for i in range(nn): 
                line = ""
                if elem_type in ["bar221i","tri231i","quad244i","quad289i"]:
                    line = f"{i+1:^10d}{reacs[j]:^10.2f}{reacs[j+1]:^10.2f}"
                    j = j+2    
                elif elem_type in ["tet344i","cube388i"]:
                    line = f"{i+1:^10d}{reacs[i]:^10.2f}{reacs[i+1]:^10.2f}{reacs[i+2]:^10.2f}"
                    j = j+3
                elif elem_type in ["beam"]:
                    line = (f"{i+1:^10d}{reacs[j]:^10.2f}{reacs[j+1]:^10.2f}{reacs[j+2]:^10.2f}"
                                     f"{reacs[j+3]:^10.2f}{reacs[j+4]:^10.2f}{reacs[j+5]:^10.2f}")
                    j = j+6

                ofl.write(line+"\n")

            ofl.write("\n\n\n")

def echo_input_file(ofl, ifile):
    ofl.write("-----------ECHO INPUT FILE-------------\n")
    with open(ifile,'r') as ifl:
        ofl.write(ifl.read())
