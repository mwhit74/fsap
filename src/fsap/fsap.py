from fsap.preprocessor import pre
from fsap.fe_processor import fep
from fsap.postprocessor import post

def fem_fi(ifile):

    """FEM by File Input

    API to call FEM program with text file input"""

    inp = pre(ifile)

    print(inp)
    
    K,Kev,F,u = fep(inp)
    
    strains,stresses,reacs = post(inp,K,Kev,F,u,ifile)

    return K,Kev,F,u,strains,stresses,reacs

def fem_di(inp):

    """FEM by Direct Input

    API to call FEM program with direct data structure input"""

    K,Kev,F,u = fep(inp)
    
    strains,stresses,reacs = post(inp,K,Kev,F,u,ifile)

    return K,Kev,F,u,strains,stresses,reacs


if __name__ == "__main__":

    ifile = argv[1]

    K,Kev,F,u,strains,stresses,reacs = fem_fi(ifile)
