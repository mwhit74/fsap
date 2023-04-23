from preprocessor import pre
from fe_processor import fep
from postprocessor import post

def gen_input_file():

    pass

def gen_dir_input():

    pass

if __name__ == " __main__":

    inp = gen_dir_input() 

    K,F,q = fep(inp)

    post(inp,K,F,q,ifile,ofile)
