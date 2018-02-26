import numpy as np

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
        for j in range(2, num_scj+1):
            if sup(i,j) == 1:
                num_reacs += 1
    ndof = num_scj*num_jt - num_reacs
    return ndof

def scv(sup, num_scj, num_jt, num_dof):
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
        str_cv (numpy array): 1D array (vector) of the numbered DOFs for
                              the entire structure to the joints
    """

    str_cv = np.empty(num_jt*num_scj)
    j = 0 #numbers for unrestrained DOFs
    k = num_dof #number for restrained DOFs
    for i in range(num_jt):
        #if a joint is in the support matrix
        #the numbering of structure coordinates changes
        #otherwise it is handled in a straight forward manner
        if i in sup[0]:
            for x in range(num_scj):
                #location of structure coordinate in vector
                #it works but it does not provide much flexibility
                y = (i - 1)*num_scj + x
                #if a joint has a support the numbers for restrained DOFs
                #are used for the restrained DOFs and otherwise the numbers
                #for unrestrained DOFs are used
                if sup[i][x] == 1:
                    str_cv[y] = k
                    k += 1
                else:
                    str_cv[y] = j
                    j += 1
        else:
            for x in range(num_scj):
                y = (i - 1)*num_scj + x
                str_cv[y] = j
                j += 1

     return str_cv

def assemble_stiffness(num_dof, num_scj, elem):
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
                    structure
        num_scj (int): number of structure coordinates per joint; number of
                       degrees of freedom per joint
        elem (list): user input member connectivity, material property
                     assignment, and section property assignment

    Returns:
        s (numpy array): assembled structure stiffness matrix
    """

    s = np.empty(num_dof, num_dof)
    gk = np.zeros(2*num_scj, 2*num_scj)
    for im in range(len(elem)):
        jb = elem[im][0]
        je = elem[im][1]
        sp_id = elem[im][2]
        a = sect[sp_id]
        matl_id = elem[im][3]
        e = matl[matl_id]
        xb = jt[im][0]
        yb = jt[im][1]
        xe = jt[im][0]
        ye = jt[im][1]
        bl = math.sqrt(math.exp(xe-xb,2) + math.exp(ye-yb,2))
        c = (xe-xb)/bl
        s = (ye-yb/bl

        mstiffg(e,a,bl,c,s,gk)
        stores()

def mstiffg(e,a,bl,c,s,gk):
    """Member global stiffness matrix."""

    eal = e*a/bl
    c2 = math.pow(c,2)
    s2 = math.pow(s,2)
    sc = s*c

    z1 = eal*c2
    z2 = eal*s2
    z3 = eal*sc

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

def stores(jb, je, num_scj, num_dof, scv, gk, s):
    for i in range(2*num_scj):
        if i <= num_scj:
            y = (jb - 1)*num_scj + 1
        else:
            y = (je - 1)*num_scj + (i - num_scj)

