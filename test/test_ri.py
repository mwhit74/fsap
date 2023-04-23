from fsap.read_input import read_input

ifile = "/home/m/Documents/python/fsap/test/bar221i_ex1b.txt"

(elem_type, const, nn, nelx, nm, nx, ns, nc, ecm, mp, sp, emp, sm, 
    njlc, jlcm, ntlc, tlcm, ndlc, dlcm, ntrlc, trlcm, blcm) = read_input(ifile)
