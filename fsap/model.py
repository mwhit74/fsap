import numpy as np
import math

def num_dof(sup, num_scj, num_jt):
    """Calculate the number of degrees of freedom of the structure.

    Loops through the user input support matrix counting the number of
    restrainted DOFs, i.e. the 1's in the 1 and 2 columns, this is the
    number of restrained coordinates. 

    The total number of DOFs is comupted as the number of joints in the
    structure which is the number of entries in the user input joint
    matrix, or coordinate matrix, times the number of DOFs per joint,
    or the number of structure coordinate joints.

    The number of unrestrained, or free,  DOFs is then calculated as the
    total number of DOFs minus the number of restrained DOFs. 
    
    Args:
        sup (list): 2D array of joint number, x support condition, y
                          support condition
        num_scj (int): number of structure coordinates per joint; number of
                       degrees of freedom per joint
        num_jt (int): number of joints in the structure

    Returns:
        ndof (int): number of unrestrained degrees of freedom of the
                    structure

    Notes:
        num_reacs (int): number of support reactions of the structure
        num_scj: truss=2, 2D frame=3, 3D frame=6
    """
    num_reacs = 0
    for i in range(len(sup)):
        for j in range(1, num_scj+1):
            if sup[i][j] == 1:
                num_reacs += 1
    ndof = num_scj*num_jt - num_reacs
    return ndof

def str_coord_vector(sup, num_scj, num_jt, ndof):
    """Generate structure coordinate vector.

    The structure coordinate vector is generated based on the number of
    restrained and unstrained DOFs. The purpose is to assign the DOF numbers
    to the correct joints and keep track of the assignment. 

    The restrained and unrestrained DOFs can also be though of as structure
    coordinates.

    The number of rows is equal to the number of structure coordinates per
    joint times the number of joints in the structure. 

    The structure coordinates relate the member stiffness coefficient matrix
    (in global coordinates) to their corresponding location in the structure
    stiffness coefficient matrix. For members meeting at a joint the
    stiffnesses of each member will be added together in structure stiffness
    matrix for each unrestrained DOF at that joint. 
    
    Args:
        sup (list): 2D array of joint number, x support condition, y support
                    condition
        num_scj (int): number of structure coordinates per joint; number of
                       degrees of freedom per joint
        num_jt (int): number of joints in the structure
        ndof (int): number of unrestrained degrees of freedom of the
                    structure

    Returns:
        scv (numpy array): 1D array (vector) of the numbered DOFs for
                              the entire structure to the joints
    """
    scv = np.empty(num_jt*num_scj, dtype='int64')
    j = 0 #numbers for unrestrained DOFs
    k = ndof #number for restrained DOFs
    # i is the assigned joint number
    for i in range(1,num_jt+1):
        #if a joint is in the support matrix
        #the numbering of structure coordinates changes
        #otherwise it is handled in a straight forward manner
        is_support = False
        for z in range(len(sup)):
            if i == sup[z][0]:
                is_support = True

                for x in range(1,num_scj+1):
                    #location of structure coordinate in vector
                    #it works but it does not provide much flexibility
                    str_coord_index = (i-1)*num_scj + x - 1
                    #if a joint has a support the numbers for restrained DOFs
                    #are used for the restrained DOFs and otherwise the numbers
                    #for unrestrained DOFs are used
                    if sup[z][x] == 1:
                        k += 1
                        scv[str_coord_index] = k
                    else:
                        j += 1
                        scv[str_coord_index] = j

        if not is_support:
            for x in range(1,num_scj+1):
                str_coord_index = (i-1)*num_scj + x - 1
                j += 1
                scv[str_coord_index] = j

    return scv

def member_properties(im, elem, matl, sect, jt):

    jb = elem[im][0]
    je = elem[im][1]
    matl_id = elem[im][2]
    e = matl[matl_id-1] #adjust index location for list
    sp_id = elem[im][3]
    a = sect[sp_id-1] #adjust index location for list
    xb = jt[jb-1][0]
    yb = jt[jb-1][1]
    xe = jt[je-1][0]
    ye = jt[je-1][1]
    bl = math.sqrt(math.pow((xe-xb),2) + math.pow((ye-yb),2))
    co = (xe-xb)/bl
    si = (ye-yb)/bl

    return jb, je, e, a, xb, yb, xe, ye, bl, co, si


def assemble_structure_stiffness_matrix(ndof, num_scj, scv, jt, 
                                        matl, sect, elem):
    """Assemble structure stiffness matrix.
   
    Gets the member connectivity data, member material property assignment,
    and member section property assignment from the 'elem' matrix. 

    Looks up the coordinates of the member begin and end joints. Calculates
    the length of the member, the sin and cos of the angle wrt the
    horizontal. 
    
    The length and sin/cos information will be used along with
    the material and section property data to generate the member stiffness
    matrix in global coordinates. The mstiffg function handles the assembly
    of of the member stiffness matrix in global coordinates.

    



    Args:
        ndof (int): number of unrestrained degrees of freedom of the
                    struture
        num_scj (int): number of structure coordinates per joint; number of
                       degrees of freedom per joint
        elem (list): user input member connectivity, material property
                     assignment, and section property assignment

    Returns:
        s (numpy array): assembled structure stiffness matrix
    """

    s = np.zeros([ndof, ndof]) #structure stiffness matrix
    gk = np.zeros([2*num_scj, 2*num_scj]) #member stiffness matrix in global
    for im in range(len(elem)):
        (jb, je, e, a, xb, yb, 
        xe, ye, bl, co, si) = member_properties(im, elem, matl, sect, jt)
        member_global_stiffness_matrix(e,a,bl,co,si,gk)
        print e
        print a
        print bl
        print co
        print si
        print gk
        store_structure_stiffness_matrix(jb, je, num_scj, ndof, scv, gk, s)

    return s

def member_global_stiffness_matrix(e,a,bl,co,si,gk):
    """Assemble member global stiffness matrix.
    
    Args:
        e (float): modulus of elasticity
        a (float): cross-section area
        bl (float): member length
        co (float): cosine of angle between member and horizontal; + ccw from
                   x-axis
        si (float): sine of angle between member and horizontal; + ccw from
                    x-axis
        gk (numpy array): member global stiffness matrix; initialized as 0
                          matrix with dimensions of member stiffness matrix

    Returns:
        None

    Notes:
        Populates member gloabl stiffness matrix by over-write zero entries.
    """

    eal = e*a/bl
    co2 = math.pow(co,2)
    si2 = math.pow(si,2)
    sico = si*co

    z1 = eal*co2
    z2 = eal*si2
    z3 = eal*sico

    gk[0][0] = z1
    gk[1][0] = z3
    gk[2][0] = -z1
    gk[3][0] = -z3
    
    gk[0][1] = z3
    gk[1][1] = z2
    gk[2][1] = -z3
    gk[3][1] = -z2

    gk[0][2] = -z1
    gk[1][2] = -z3
    gk[2][2] = z1
    gk[3][2] = z3

    gk[0][3] = -z3
    gk[1][3] = -z2
    gk[2][3] = z3
    gk[3][3] = z2

def store_structure_stiffness_matrix(jb, je, num_scj, ndof, scv, gk, s):
    """Adds elements of member global stiffness to structure stiffness matrix.

    Step down rows first then across columns. The rows represent the forces in
    the directions of the DOFs. The columns represent the displacements in the
    directions of the DOFs. 

    When the row number is less than the number of coordinates per joint it is
    representing a DOF at the beginning of a member. When the row number is
    greater a DOF at then end of a member. The index location in the structure
    coordinate vector is calculated based on the row number. The code number 
    which represents the DOF of the structure is located with the calculated 
    index location. If the code number is less than the number of unrestrained
    DOFs then a similar process is carried out for the column number. If it
    found that code number for the column is also less than the number of
    restrained DOFs then the member global stiffness coefficient at the local of
    the row and column is added to the structure stiffness coefficient located 
    at the associated code number location. 

    Args:
        jb (int): begin joint id
        je (int): end joint id
        num_scj (int): number of structure coordinates per joint; number of
                       degrees of freedom per joint
        ndof (int): number of unrestrained degrees of freedom of the
                       structure
        scv (numpy array): 1D array (vector) of the numbered DOFs for
                              the entire structure to the joints
        gk (numpy array): 2D array representing global member stiffness matrix
        s (numpy array): 2D array representing assembled structure stiffness 
                         matrix; initialized as empty array with dimensions of
                         the structure stiffness matrix
    Returns:
        None

    Notes:
        Populates structure stiffness matrix
    """
    #stepping down rows
    #i is the row
    for member_coord_row in range(1,2*num_scj+1):
        #if i <= num_scj then i member coordinate for begin joint
        #else i is member coordinate for end joint 
        if member_coord_row <= num_scj:
            str_coord_row_index = (jb - 1)*num_scj + member_coord_row - 1
        else:
            #(i - num_scj) is accounting for the offset
            str_coord_row_index = ((je - 1)*num_scj + 
                                    (member_coord_row - num_scj - 1))
        str_coord_row = scv[str_coord_row_index] 
        #if n1 is an unrestrained DOF
        if str_coord_row <= ndof:
            #stepping across rows
            #j is the column
            for member_coord_col in range(1,2*num_scj+1):
                if member_coord_col <= num_scj:
                    str_coord_col_index = ((jb - 1)*num_scj + 
                                            member_coord_col - 1)
                else:
                    str_coord_col_index = ((je - 1)*num_scj + 
                                            (member_coord_col - num_scj - 1))
                str_coord_col = scv[str_coord_col_index]
                #if n2 is an unrestrained DOF
                if str_coord_col <= ndof:
                    s[str_coord_row-1][str_coord_col-1] = (s[str_coord_row-1][str_coord_col-1] +
                                                      gk[member_coord_row-1][member_coord_col-1])


def joint_load_vector(ndof, num_scj, scv, load):
    """Assemble the joint load vector.


    """
    p = np.zeros(ndof)
    for x in range(len(load)):
        jt = load[x][0]
        str_coord_index = (jt - 1)*num_scj - 1
        for load_index in range(num_scj):
            str_coord_index = str_coord_index + 1
            str_coord = scv[str_coord_index]
            if str_coord <= ndof:
                p[str_coord-1] = p[str_coord-1] + load[x][load_index+1]
    
    return p


def member_forces_disps_reacs(num_dof, num_scj, scv,
                              elem, matl, sect, jt, gdisp):
    """Calculate the member forces and displacements, and support reactions.


    """
    mlsm = np.zeros([2*num_scj,2*num_scj]) #member local stiffness matrix
    mtm = np.zeros([2*num_scj,2*num_scj]) #member transform matrix
    mged = np.zeros(2*num_scj) #member global end disps
    mled = np.zeros(2*num_scj) #member local end disps
    mlef = np.zeros(2*num_scj) #member local end forces
    mgef = np.zeros(2*num_scj) #member global end forces
    reac = np.zeros(2*num_scj) #support reactions

    for im in range(len(elem)):
        (jb, je, e, a, xb, yb, 
        xe, ye, bl, co, si) = member_properties(im, elem, matl, sect, jt)
        print "\nmember: " + str(im)
        member_global_disp(jb, je, num_scj, num_dof, scv, gdisp, mged)
        print mged
        member_transform_matrix(co, si, num_scj, mtm)
        print mtm
        member_local_disp(num_scj, mged, mtm, mled)
        print mled
        member_local_stiffness_matrix(e, a, bl, num_scj, mlsm)
        print mlsm
        #member_local_forc()
        #member_global_forc()
        #sup_reac()


def member_global_disp(jb, je, num_scj, num_dof, scv, gdisp, mged):
    """Member global displacement vector.

    For each member, transfers the displacment in global coordinates from 
    the structure coordinates to the local member in global coordinates.

    Args:
        jb
        je
        num_scj
        num_dof
        scv
        gdisp(numpy array): displacement of the structure coordinates in global
                            coordinates
        mged(numpy array): local member end displacements in global coordinates

    Returns:
        None

    Notes:
        Populates the mged array with values.
    """
    #look at beginning joint of each member
    str_coord_index = (jb - 1)*num_scj - 1
    for i in range(num_scj):
        str_coord_index = str_coord_index + 1
        str_coord = scv[str_coord_index]
        if str_coord <= num_dof:
            mged[i] = gdisp[str_coord-1]

    #look at end joint of each member
    str_coord_index = (je - 1)*num_scj - 1
    for i in range(num_scj,2*num_scj):
        str_coord_index = str_coord_index + 1
        str_coord = scv[str_coord_index]
        if str_coord <= num_dof:
            mged[i] = gdisp[str_coord-1]


def member_transform_matrix(co, si, num_scj, mtm):
    """Assemble member transformation matrix.
    
    Stores the values of the direction sines and cosines in the appropriate
    places in the matrix.

    Assembles the member local to global transformation matrix

    Args:
        co(float): cosine of the angle between horizontal and member axis;
                   global x-axis positive to right
        si(float): sine of the angle between horizontal and member axis; global
                   x-axis positive to right
        num_scj
        mtm(numpy array): 2D arry representation of the transformation matrix

    Returns:
        None

    Notes:
        Populates the mtm matrix.
    """
    mtm[0][0] = co
    mtm[0][1] = si
    mtm[1][0] = -si
    mtm[1][1] = co
    mtm[2][2] = co
    mtm[2][3] = si
    mtm[3][2] = -si
    mtm[3][3] = co


def member_local_disp(num_scj, mged, mtm, mled):
    """Member displacements in local coodinates.

    """
    for i in range(2*num_scj):
        for j in range(2*num_scj):
            mled[i] = mled[i] + mtm[i][j]*mged[j]


def member_local_stiffness_matrix(e, a, bl, num_scj, mlsm):
    """Assemble member stiffness matrix in local coordinates.

    """
    for i in range(2*num_scj):
        for j in range(2*num_scj):
            z = e*a/bl
            mlsm[0][0] = z
            mlsm[0][2] = -z
            mlsm[2][0] = -z
            mlsm[2][2] = z



def member_local_forc():
    pass


def member_global_forc():
    pass


def sup_reac():
    pass


