from fsap.read_input import read_input

def pre(ifile):

    (elem_type, const, nn, nelx, nm, nx, ns, nc, ecm, mp, sp, emp, epa, sm, 
        njlc, jlcm, ntlc, tlcm, ndlc, dlcm, ntrlc, trlcm, blcm, wob, ofile) = read_input(ifile)

    #print(read_input(ifile))

    return (elem_type, const, nn, nelx, nm, nx, ns, nc, ecm, mp, sp, 
            emp, epa, sm, njlc, jlcm, ntlc, tlcm, ndlc, dlcm, ntrlc, trlcm, blcm, wob, ofile)

def mesh():

    pass


