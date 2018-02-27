"""Run, Solve, Execute!"""

import ri
import model
import naf.linalg.lu 

input_file = r''

jt, sup, matl, sect, elem, load = ri.read_input(input_file)

num_scj = 2 #truss
num_jt = len(jt)

ndof = num_dof(sup, num_scj, num_jt)
scv = str_coord_vector(sup, num_scj, num_jt, ndof)
sk = assemble_stiffness(ndof, num_scj, elem)
p = joint_load_vector(ndof, num_scj, scv, load)

LU, ov = lu.lu_decomp(sk)
x = lu.lu_solve(LU, ov, p)


