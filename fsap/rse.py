"""Run, Solve, Execute!"""

import ri
import model
import output
import naf.linalg.lu as lu

input_file = ('test/sample_input_file.txt')

jt, sup, matl, sect, elem, load = ri.read_input(input_file)

num_scj = 2 #truss
num_jt = len(jt)

ndof = model.num_dof(sup, num_scj, num_jt)
scv = model.str_coord_vector(sup, num_scj, num_jt, ndof)
sk = model.assemble_structure_stiffness_matrix(ndof, num_scj, scv, jt, matl, sect, elem)
p = model.joint_load_vector(ndof, num_scj, scv, load)

LU, ov = lu.lu_decomp(sk)
gdisp = lu.lu_solve(LU, ov, p)

(mgedl, mtml, mledl, mlsml, 
mlefl, mgefl, reac) = model.member_forces_disps_reacs(ndof, num_scj, scv,
                                                      elem, matl, sect, jt, gdisp)
output.echo_input(jt, sup, matl, sect, elem, load)
output.echo_output(gdisp, mlefl, reac, scv, jt, num_scj, ndof)
