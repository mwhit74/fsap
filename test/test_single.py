from fsap.preprocessor import pre
from fsap.fe_processor import fep
from fsap.postprocessor import post
import pathlib

ifile = "/home/m/Documents/python/fsap/test/bar221i_ex1_opt.txt"

ipath = pathlib.Path(ifile)
ofile = pathlib.Path.joinpath(ipath.parent,ipath.stem+"_out.txt")

(elem_type, const, nn, nelx, nm, nx, ns, nc, ecm, mp, sp, 
            emp, sm, njlc, jlcm, ntlc, tlcm, ndlc, dlcm, ntrlc, trlcm, blcm) = pre(ifile)

#print(elem_type, const, nn, nelx, nm, nx, ns, nc, ecm, mp, sp,
#          emp, sm, njlc, jlcm, ntlc, tlcm, ndlc, dlcm, ntrlc, trlcm, blcm)

inp = pre(ifile)

K,F,q = fep(inp)

#print(K,F,q)
post(inp,K,F,q,ifile,ofile)


