from fsap.preprocessor import pre
from fsap.fe_processor import fep
from fsap.postprocessor import post
from fsap.fsap import fem_fi,fem_di
import pathlib

ifiles = ["/home/m/Documents/python/fsap/test/bar221i_ex1.txt",
          "/home/m/Documents/python/fsap/test/bar221i_ex2.txt",
          "/home/m/Documents/python/fsap/test/tri231i_ex1.txt"]


for ifile in ifiles:

    ipath = pathlib.Path(ifile)

    K,Kev,F,u,strains,stresses,reacs = fem_fi(ifile)

    #(elem_type, const, nn, nelx, nm, nx, ns, nc, ecm, mp, sp, 
    #            emp, epa, sm, njlc, jlcm, ntlc, tlcm, ndlc, dlcm, ntrlc, trlcm, blcm, wob, ofile) = pre(ifile)
    
    #print(elem_type, const, nn, nelx, nm, nx, ns, nc, ecm, mp, sp,
    #          emp, epa, sm, njlc, jlcm, ntlc, tlcm, ndlc, dlcm, ntrlc, trlcm, blcm, wob, ofile)
    
    #inp = pre(ifile)
    #
    #K,Kev,F,u = fep(inp)
    #
    #strains, stress, reacs = post(inp,K,Kev,F,u,ifile)


