import math
import numpy as np
from numpy import matmul, transpose, linalg, outer
                 
def fep(inp):

    """Finite Element Processor

    Params:
        inp (tuple): input data

    R (for optimization)eturns:
        K,F,u (tuple)
            K (ndarray): global stiffness matrix w/ BCs
            F (ndarray): global load vector w/ BCs
            u (ndarray): nodal displacements

    Defs:
        elem_type (str): element type
        const (int): constitutive relationship id
        nn (int): number of global nodes
        nelx (int): number of global elements
        nm (int): number of materials
        nx (int): number of cross-sections
        ns (int): number of supports
        nc (ndarray): nodal coordinate matrix
        ecm (ndarray): element connection matrix
        mp (ndarray): material property matrix
        sp (ndarray): section properties matrix
        epm (ndarray): assigned element property matrix
        epa (ndarray): design variables for optimization
        sm (ndarray): support matrix
        njlc (int): number of joint load cases
        jlcm (ndarray): joint load case matrix
        ntlc (int): number of temperature load cases
        tlcm (ndarray): temperature load case matrix
        ntrlc (int): number of something load cases
        trlcm (ndarray): something load case matrix
        wob (bool): flag to write output file
        ofile (str): path to output file
    """

    (elem_type, const, nn, nelx, nm, nx, ns, nc, ecm, mp, sp, 
            epm, epa, sm, njlc, jlcm, ntlc, tlcm, ndlc, dlcm, 
                            ntrlcm, trlcm, blcm, wob, ofile) = inp

    #assemble_type = 1
    #Kev - vector of element stiffness matrices in global (for optimization)
    K,Kev = assemble_global_stiffness_matrix(elem_type, const, nn, nelx, nc,
            ecm, mp, sp, epm, epa)

    F = assemble_load_vector(elem_type,nn,jlcm,tlcm,trlcm,blcm)

    Kmod = K.copy()
    Fmod = F.copy()
    #bc_type = 1
    Kmod,Fmod = apply_bcs(elem_type,Kmod,Fmod,nn,ns,sm,ndlc,dlcm)

    soln_type = 1
    u = solve_equs(Kmod,Fmod,soln_type=1)

    return K,Kev,F,u


def solve_equs(K,F,soln_type=1):

    """Solve Linear System of Equations

    Params:
        K (ndarray): global stiffness matrix w/ BCs
        F (ndarray): global load vector w/ BCs
        soln_type (int,optional): solution method type, default=1

    Returns:
        u (ndarray): nodal displacements

    Notes:
        1. soln_type options:
            1) general solver
            2)

    """

    if soln_type == 1:
        u = linalg.solve(K,F)

    return u

def assemble_global_stiffness_matrix(elem_type, const, nn, nelx, nc, ecm,
                                         mp, sp, epm, epa, assemble_type=1):

    """Assemble Global Stiffness Matrix

    Params:
        elem_type (str): unique element id
        const (int): constitutive relationship id
        nn (int): number of global nodes
        nelx (int): number of global elements
        nc (ndarray): nodal coordinate matrix
        ecm (ndarray): element connection matrix
        mp (ndarray): material property matrix
        sp (ndarray): section properties
        epm (ndarray): assigned element property matrix
        epa (ndarray): design variables for optimization
        assemble_type (int,optional): method used to assemble global 
                                      stiffness matrix, default=1

    Returns:
        K (ndarray): assembled global stiffness matrix

    Notes:
        1. assemble_type options:
            1) full matrix assembly
            2) reduced bandwidth
            3) skyline

    """

    ndofn = elementdata[elem_type][0]
    nne = elementdata[elem_type][1] # number of nodes per element
    nedof = ndofn*nne # total number of dofs per element

    Kev = np.empty((nelx,nedof,nedof))
    
    if assemble_type == 1:
        ngdof = ndofn*nn #number of global DOFs
        K = np.zeros([ngdof, ngdof])
        for i in range(nelx):

            nodes = ecm[i,1:]
            elem_fnc = elementdata[elem_type][3]
            dofs = calc_dofs(elem_type,nodes)

            Ke = assemble_local_stiffness_matrix(i,elem_type,const,
                                                 nc,ecm,mp,sp,epm,epa,elem_fnc)

            Kev[i,:,:] = Ke

            for j in range(dofs.shape[0]):
                f = dofs[j]
                for k in range(dofs.shape[0]):
                    g = dofs[k]
                    K[f,g] = K[f,g] + Ke[j,k]

    return K,Kev

def calc_dofs(elem_type, nodes):

    """Calculate global DOFs for element

    Params:
        elem_type (str): element type
        nodes (ndarray): global element node numbers

    Returns:
        dofs (ndarray): global dof numbers
    """
    ndofn = elementdata[elem_type][0]
    nne = elementdata[elem_type][1]

    m = 0
    dofs = np.empty((nne*ndofn),dtype='int')

    for i in range(nne):
        k = ndofn - 1
        n = nodes[i]
        for j in range(ndofn):
            dofs[m] = ndofn*n - k - 1
            k=k-1
            m=m+1

    return dofs
    
def assemble_local_stiffness_matrix(i,elem_type,const,nc,ecm,mp,sp,epm,epa,elem_fnc):

    """Assemble local stiffness matrix

    Params:
        elem_type (str): unique element id
        const (int): constitutive relationship id
        nc (ndarray): nodal coordinate matrix
        ecm (ndarray): element connection matrix
        mp (ndarray): material property matrix
        sp (ndarray): section properties
        epm (ndarray): assigned element property matrix
        epa (ndarray): design variables for optimization
        elem_fnc (str): element stiffness matrix function name

    Returns:
        Ke (ndarray): local element stiffness matrix in global coords

    """
    ndofn = elementdata[elem_type][0] # number of dofs per node for element
    nne = elementdata[elem_type][1] # number of nodes per element
    nip = elementdata[elem_type][2] # number of integration points
    qt_dict = elementdata[elem_type][4] # gauss quad dict for element
    qt_id = elementdata[elem_type][5] # gauss quad value id

    qt = qt_dict[qt_id]
    ips = qt[0] # gauss quad locations
    wts = qt[1] # gauss quad weights

    nedof = ndofn*nne # total number of dofs per element

    Ke = np.zeros((nedof,nedof))

    for j in range(nip):
        
        wt = wts[j] 

        ip = ips[j]

        S,D,B = elem_fnc(i,const,nc,ecm,mp,sp,epm,epa,ip)

        #bar elements
        if D.shape[0] == 1:

            Ke = Ke + wt*D[0]*outer((B),B)*S

        #everything else
        else:
 
            Ke = Ke + wt*matmul(transpose(B),matmul(D,B))*S

    return Ke
        

def bar221i(i,const,nc,ecm,mp,sp,epm,epa,ip):

    """2D Linear Isoparametric Bar Element data

    Params:
        i (int): element number
        const (int): constitutive relationship id
        nc (ndarray): nodal coordinate matrix
        ecm (ndarray): element connectivity matrix
        mp (ndarray): material property matrix
        sp (ndarray): section property matrix
        epm (ndarray): element property assignment matrix
        epa (ndarray): design variables for optimization
        ip (ndarray): integration point information, not used

    Returns:
        Ke (ndarray): element stiffness matrix in global coordinates

    """

    n1 =ecm[i,1]-1 #minus 1 converts to index position
    n2 =ecm[i,2]-1 #minus 1 converts to index position

    x1 = nc[n1,1]
    y1 = nc[n1,2]

    x2 = nc[n2,1]
    y2 = nc[n2,2]

    le = math.sqrt((x2-x1)**2+(y2-y1)**2)
    l = (x2-x1)/le
    m = (y2-y1)/le

    DJ = le/2

    B = 1/le*np.array([-1,1])

    L = np.array([[l,m,0,0],[0,0,l,m]])

    B = matmul(B,L)

    sp_num = epm[i,1]-1 #minus 1 to offset index from user input
    dv = epa[i,1] #design variable value
    Area = dv*sp[sp_num,1] #multiply area by design variable for optimization code
    S = Area*DJ*2 #the 2 comes from the Guass quadrature integration

    D = calc_D(const,i,mp,epm,epa) #returns E, modulus of elasticity

    return S,D,B

def bar321i(i,const,nc,ecm,mp,sp,epm,epa,ip):

    """3D Linear Isoparametric Bar Element data

    Params:
        i (int): element number
        const (int): constitutive relationship id
        nc (ndarray): nodal coordinate matrix
        ecm (ndarray): element connectivity matrix
        mp (ndarray): material property matrix
        sp (ndarray): section property matrix
        epm (ndarray): element property assignment matrix
        epa (ndarray): design variables for optimization
        ip (ndarray): integration point information, not used

    Returns:
        Ke (ndarray): element stiffness matrix in global coordinates

    """

    n1 = ecm[i,1]-1 #minus 1 converts to index position
    n2 = ecm[i,2]-1 #minus 1 converts to index position
    n3 = ecm[i,3]-1 #minus 1 converts to index position

    x1 = nc[n1,1]
    y1 = nc[n1,2]
    z1 = nc[n1,3]

    x2 = nc[n2,1]
    y2 = nc[n2,2]
    z2 = nc[n2,3]

    le = math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    l = (x2-x1)/le
    m = (y2-y1)/le
    n = (z2-z1)/le

    DJ = le/2

    B = 1/le*np.array([-1,1])

    L = np.array([[l,m,n,0,0,0],[0,0,0,l,m,n]])

    B = matmul(B,L)

    sp_num = epm[i,1]-1 #minus 1 to offset index from user input, 
    dv = epa[i,1] #design variable value
    Area = dv*sp[sp_num,1] #multiply area by design variable for optimization code
    S = Area*DJ*2 #the 2 comes from the Guass quadrature integration

    D = calc_D(const,i,mp,epm) #returns E, modulus of elasticity

    return S,D,E

def tri231i(i,const,nc,ecm,mp,sp,epm,epa,ip):

    """2D 3-node Isoparametric Linear Triangular Element

    Params:
        i (int): element number
        const (int): constitutive relationship id
        nc (ndarray): nodal coordinate matrix
        ecm (ndarray): element connectivity matrix
        mp (ndarray): material property matrix
        sp (ndarray): section property matrix
        epm (ndarray): element property assignment matrix
        epa (ndarray): design variables for optimization
        ip (ndarray): integration point information, not used

    Returns:
        Ke (ndarray): element stiffness matrix in global coordinates

    """

    n1 = ecm[i,1]-1 #minus 1 converts to index position
    n2 = ecm[i,2]-1 #minus 1 converts to index position
    n3 = ecm[i,3]-1 #minus 1 converts to index position

    x1 = nc[n1,1]
    y1 = nc[n1,2]

    x2 = nc[n2,1]
    y2 = nc[n2,2]
    
    x3 = nc[n3,1]
    y3 = nc[n3,2]

    DJ = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)

    B = 1/DJ*np.array([[y2-y3,0,y3-y1,0,y1-y2,0],
                       [0,x3-x2,0,x1-x3,0,x2-x1],
                       [x3-x2,y2-y3,x1-x3,y3-y1,x2-x1,y1-y2]])

    sp_num = epm[i,1]-1 #minus 1 to offset index from user input
    t = sp[sp_num,1]
    Area = DJ/2
    S = Area*t

    D = calc_D(const,i,mp,epm)

    return S,D,B

def tri263i(i,const,nc,ecm,mp,sp,epm,epa,ip):

    """2D 6-node Isoparametric Linear Triangular Element

    Params:
        i (int): element number
        const (int): constitutive relationship id
        nc (ndarray): nodal coordinate matrix
        ecm (ndarray): element connectivity matrix
        mp (ndarray): material property matrix
        sp (ndarray): section property matrix
        epm (ndarray): element property assignment matrix
        epa (ndarray): design variables for optimization
        ip (ndarray): integration point information

    Returns:
        Ke (ndarray): element stiffness matrix in global coordinates

    """
    n1 = ecm[i,1]-1 #minus 1 converts to index position
    n2 = ecm[i,3]-1 #minus 1 converts to index position
    n3 = ecm[i,5]-1 #minus 1 converts to index position
    n4 = ecm[i,2]-1 #minus 1 converts to index position
    n5 = ecm[i,4]-1 #minus 1 converts to index position
    n6 = ecm[i,6]-1 #minus 1 converts to index position

    x1 = nc[n1,1]
    y1 = nc[n1,2]

    x2 = nc[n2,1]
    y2 = nc[n2,2]
    
    x3 = nc[n3,1]
    y3 = nc[n3,2]

    x4 = nc[n4,1]
    y4 = nc[n4,2]

    x5 = nc[n5,1]
    y5 = nc[n5,2]

    x6 = nc[n6,1]
    y6 = nc[n6,2]
    
    x = np.array([x1,x2,x3,x4,x5,x6])
    y = np.array([y1,y2,y3,y4,y5,y6])

    xi1 = ip[0]
    xi2 = ip[1]
    xi3 = 1-xi1-xi2
    dNdxi = np.array([4*xi1-1,0,-4*xi3+1,4*xi2,-4*xi2+4*(xi3-xi1)])
    dNdeta = np.array([0,4*xi2-1,-4*xi3+1,4*xi2,4*(xi3-xi2),4*xi3])

    S,D,B = calc_SDB_2D(x,y,dNdxi,dNdeta,const,i,mp,epm)
    
    return S,D,B

def quad244i(i,const,nc,ecm,mp,sp,epm,epa,ip):
    """2D 4-node Isoparametric Linear Quad Element

    Params:
        i (int): element number
        const (int): constitutive relationship id
        nc (ndarray): nodal coordinate matrix
        ecm (ndarray): element connectivity matrix
        mp (ndarray): material property matrix
        sp (ndarray): section property matrix
        epm (ndarray): element property assignment matrix
        epa (ndarray): design variables for optimization
        ip (ndarray): integration point information

    Returns:
        Ke (ndarray): element stiffness matrix in global coordinates

    """
    n1 = ecm[i,0]-1 #minus 1 converts to index position
    n2 = ecm[i,1]-1 #minus 1 converts to index position
    n3 = ecm[i,2]-1 #minus 1 converts to index position
    n4 = ecm[i,3]-1 #minus 1 converts to index position

    x1 = nc[n1,1]
    y1 = nc[n1,2]

    x2 = nc[n2,1]
    y2 = nc[n2,2]
    
    x3 = nc[n3,1]
    y3 = nc[n3,2]

    x4 = nc[n4,1]
    y4 = nc[n4,2]

    x = np.array([x1,x2,x3,x4])
    y = np.array([y1,y2,y3,y4])

    xi = ip[0]
    eta = ip[1]
    dNdxi = np.array([-0.25*(1-eta),0.25*(1-eta),0.25*(1+eta),-0.25*(1+eta)])
    dNdeta = np.array([-0.25*(1-xi),-0.25*(1+xi),0.25*(1+xi),0.25,(1-xi)])

    S,D,B = calc_SDB_2D(x,y,dNdxi,dNdeta,const,i,mp,epm)
    
    return S,D,B

def quad289i():

    pass

def calc_SDB_2D(x,y,dNdxi,dNdeta,const,i,mp,epm,epa):

    """Transformed area, constitutive matrix, strain-disp matrix

    Calculates the transformed area with the Jacobian matrix, selects
    the appropriate constitutive matrix, and calculates the 
    strain-dispalcment matrix

    Params:
        x (ndarray): array of nodal x-coords
        y (ndarray): array of nodal y-coords
        dNdxi (ndarray): partial derivatives of shape functions
        dNdeta (ndarray): partical derivatives of shape functions
        const (int): constitutive relationship id
        i (int): element number
        mp (ndarray): material property matrix
        epm (ndarray): assigned element property matrix
        epa (ndarray): design variables for optimization

    Returns:
        S (ndarray): transformed area
        D (ndarray): constitutive relationship matrix
        B (ndarray): strain-displacement matrix
    """

    nsf = dNdxi.shape[0] #number of shape functions

    J = np.empty(2,2)
    J[0,0] = matmul(dNdxi,transpose(x))
    J[0,1] = matmul(dNdxi,transpose(y))
    J[1,0] = matmul(dNdeta,transpose(x))
    J[1,1] = matmul(dNdeta,transpose(y))

    DJ = J[0,0]*J[1,1]-J[1,0]*J[0,1]

    A = np.empty(3,4)
    A[0,0] = J[1,1]
    A[0,1] = -J[0,1]
    A[0,2] = 0
    A[0,3] = 0
    A[1,0] = 0
    A[1,1] = 0
    A[1,2] = -J[1,0]
    A[1,3] = J[0,0]
    A[2,0] = -J[1,0]
    A[2,1] = J[0,0]
    A[2,2] = J[1,1]
    A[2,3] = -J[0,1]

    G = np.empty(4,nsf)
    for i in range(2):
        if i == 0:
            arr = dNdxi
        else:
            arr = dNdeta
        k = 0
        for j in range(2*nsf):
            if j%2 == 0:
                v = arr[k]
                k = k+1
            else:
                v = 0
            G[i,j] = v
            G[i+1,j] = v

    B = 1/DJ*matmul(A,G)

    sp_num = epm[i,1]-1 #minus 1 to offset index from user input
    t = sp[sp_num,1]
    Area = DJ #DJ*dxi*deta
    S = t*Area

    D = calc_D(const,i,mp,epm)

    return S,D,B


def tet344i():

    pass

def cube388i():

    pass

def calc_SDB_3D(x,y,dNdxi,dNdeta,const,i,mp,epm,epa):

    pass

def calc_D(const,i,mp,epm,epa):
    
    """Calculate constitutive relationship matrix

    Params:
        const (int): constitutive relationship id
        i (int): element number
        mp (ndarray): material property matrix
        epm (ndarray): assigned element property matrix
        epa (ndarray): design variables for optimization
    
    Returns:
        D (ndarray): constitutive relatioship matrix
    """

    if const == 1:
        mat_prop_num = epm[i,1]-1 #minus 1 to offset index from user input
        D = np.array([mp[mat_prop_num,1]])

    elif const == 2:
        mat_prop_num = epm[i,1]-1 #minus 1 to offset index from user input
        E = mp[mat_prop_num,1]
        v = mp[mat_prop_num,2]

        D = plane_stress(E,v)

    elif const == 3: 
        mat_prop_num = epm[i,1]-1 #minus 1 to offset index from user input
        E = mp[mat_prop_num,1]
        v = mp[mat_prop_num,2]

        D = plane_strain(E,v)

    elif const == 4:
        mat_prop_num = epm[i,1]-1 #minus 1 to offset index from user input
        E = mp[mat_prop_num,1]
        v = mp[mat_prop_num,2]

        D = const_3d(E,v)

    return D
   
def plane_stress(E,v):

    """Plane Stress Linear Elastic Constitutive Relation

    Params:
        E (float): modulus of elasticity
        v (float): Poisson's ratio

    Returns:
        D (ndarray): Plane Stress constitutive relation

    """

    D = E/(1-v**2)*np.array([[1,v,0],
                             [v,1,0],
                             [0,0,(1-v)/2]])

    return D

def plane_strain(E,v):

    """Plane Strain Linear Elastic Constitutive Relation

    Params: E (float): modulus of elasticity
        v (float): Poisson's ratio

    Returns:
        D (ndarray): Plane Strain constitutive relation

    """
    D = E/((1+v)*(1-2*v))*np.array([[(1-v),v,0],
                                    [v,(1-v),0],
                                    [0,0,(1-2*v)/2]])

    return D

def const_3d(E,v):

    """3D Linear Elastic Constitutive Relation

    Params:
        E (float): modulus of elasticity
        v (float): Poisson's ratio

    Returns:
        D (ndarray): 3D constitutive relation

    """
    a = 1-v
    b = (1-2*v)/2

    D = E/((1+v)*(1-2*v))*np.array([[a,v,v,0,0,0],
                                    [v,a,v,0,0,0],
                                    [v,v,a,0,0,0],
                                    [0,0,0,b,0,0],
                                    [0,0,0,0,b,0],
                                    [0,0,0,0,0,b]])

    return D

def assemble_na_qt():

    """Dummy Gauss quadrature

    A dummy entry for elements that do not require integration,
    i.e. linear bar, constant strain triangle
    """

    na_qt = {0:([1.0],[1.0])}

    return na_qt

def assemble_lin_qt():

    """Linear Gauss quadrature"""

    #locations
    a = 0.5773502692
    b = 0.7745966692
    #weights
    c = 0.5555555556
    d = 0.8888888889

    lin_qt = {0:([0.0],[2.0]),
              1:([-a,a],[1.0]),
              2:([-b,0.0,b],[c,d])}

    return lin_qt

def assemble_tri_qt():

    """Triangular Gauss quadrature"""


    #locations
    a = 0.1666666666
    b = 0.6666666666
    #weights
    c = a
    
    tri_qt = {0:([(a,a),(b,a),(a,b)],[c,c,c])}

    return tri_qt

def assemble_quad_qt():

    """Quad Gauss quadrature"""

    #locations
    a = 0.5773502692
    b = 0.7745966692
    #weights
    c = 0.5555555556
    d = 0.8888888889
    #weights multiplied
    e = c**2
    f = c*d

    quad_qt = {0:([(0.0,0.0)],[2.0]), #1x1
               1:([(0.0,-a),(0.0,a)],[1.0,1.0]), #1x2
               2:([(-a,a),(a,a),(-a,-a),(a,-a)],[1.0,1.0,1.0,1.0]), #2x2
               3:([(-b,b),(0.0,b),(b,b),
                   (-b,0.0),(0.0,0.0),(b,0.0),
                   (-b,-b),(0.0,-b),(b,-b)],
                   [e,f,e,e,e,e,f,e,f])} #3x3

    return quad_qt

def tet_qt():

    pass

def cube_qt():

    pass

#Initilize Gauss quadrature dicts
na_qt = assemble_na_qt()
lin_qt = assemble_lin_qt()
tri_qt = assemble_tri_qt()
quad_qt = assemble_quad_qt()

# [element type][dim][number of nodes][number of int. pts][I-iso|S-sub]: 
#       (ndofn,nne,ip,ke_fcn,qt_dict,qt_id)
# ndofn - number of degrees of freedom per node
# nne - number of nodes per element
# ip - number of integration points
# ke_fcn - element stiffness function name
# qt_dict - quadrature dict 
# qt_id - quadrature id
# ssd - strain and stress dim
elementdata = {"bar221i":(2,2,1,bar221i,na_qt,0,1),
               "bar321i":(3,2,1,bar321i,na_qt,0,1),
               "tri231i":(2,3,1,tri231i,na_qt,0,3),
               "tri263i":(2,6,3,tri263i,tri_qt,0,3),
               "quad244i":(2,4,4,quad244i,quad_qt,2,3),
               "quad289i":(2,8,9,quad289i,quad_qt,3,3),
               "tet344i":(3,4,4,tet344i,tet_qt,0,6),
               "cube388i":(3,8,8,cube388i,cube_qt,6)}


    
def assemble_load_vector(elem_type,nn,jlcm,tlcm,trlcm,blcm):

    """Assemble Load Vector

    Params:
        nn (int): number of nodes
        jlcm (ndarray): joint load vector
        tlcm (ndarray): temperature load vector
        trlcm (ndarray): traction load vector
        blcm (ndarray): body load vector

    Returns:
        F (ndarray): load vector
    """

    ndofn = elementdata[elem_type][0]

    njlc = jlcm.shape[0]
    ntlc = tlcm.shape[0]
    ntrlc = trlcm.shape[0]
    nblc = blcm.shape[0]

    Fj = np.zeros((ndofn*nn))
    Ft = np.zeros((ndofn*nn))
    Ftr = np.zeros((ndofn*nn))
    Fb = np.zeros((ndofn*nn))
    
    if njlc != 0:
        Fj = joint_load_vector(ndofn,nn,jlcm)
    if ntlc != 0:
        Ft = temp_load_vector(ndofn,nn,blcm)
    if ntrlc != 0:
        Ftr = traction_load_vector(ndofn,nn,trlcm)
    if nblc != 0:
        Fb = body_load_vector(ndofn,nn,tlcm)
   
    F = Fj + Ft + Ftr + Fb

    return F
    

def joint_load_vector(ndofn,nn,jlcm):

    """Joint Load Vector

    Params:
        nn (int): number of nodes
        nl (int): number of loads
        lm (ndarray): load matrix

    Returns:
        F (ndarray): global load vector

    Notes:
        1. This should work for both 2D and 3D nodal loads
        2. This should work for concentrated loads and concentrated moments
    """

    F = np.zeros((ndofn*nn))
    k = ndofn 

    for i in range(jlcm.shape[0]):
        jlm = jlcm[i,1:]
        node = int(jlm[0])
        dirc = int(jlm[1])
        mag = jlm[2]

        if dirc == 11:
            dof = ndofn*node - k
            F[dof] = mag
        if dirc == 12:
            dof = ndofn*node - (k-1)
            F[dof] = mag
        if dirc == 13:
            dof = ndofn*node - (k-2)
            F[dof] = mag
        if dirc == 21:
            dof = ndofn*node - (k-3)
            F[dof] = mag
        if dirc == 22:
            dof = ndofn*node - (k-4)
            F[dof] = mag
        if dirc == 23:
            dof = ndofn*node - (k-5)
            F[dof] = mag

    return F

def temp_load_vector(ndofn,nn,tlcm):

    F = np.zeros((ndofn*nn))

    return F

def body_load_vector(ndofn,nn,blcm):

    F = np.zeros((ndofn*nn))

    return F

def traction_load_vector(ndofn,nn,trlcm):

    F = np.zeros((ndofn*nn))

    return F

def apply_bcs(elem_type,K,F,nn,ns,sm,ndlc,dlcm,bc_type=1):

    """Applies Boundary Conditions

    Params:
        elem_type (str): unique element id
        K (ndarray): global stiffness matrix w/o BCs
        F (ndarray): global load vector w/o BCs
        nn (int): number of nodes
        ns (int): number of supported nodes
        sm (ndarray): support matrix
        ndlc (int): number of support displacements
        dlcm (ndarray): support displacement matrix
        bc_type (int,optional): select method of applying BCs, default=1
        C (float,optional): penalty scaling factor, default=10E4

    Returns:
        K (ndarray): global stiffness matrix with BCs applied
        F (ndarray): global load vector with BCs applied

    Notes:
        1. bc_type options: 
            1) Penalty Method, C=10e4
            2) Elimination Method
    """

    ndofn = elementdata[elem_type][0]

    if bc_type == 1:

        Kmax = K.max()
        sc = 10E4 #stiffness scalar
        C = sc*Kmax
        K = k_penalty_bc(K,ndofn,ns,sm,C)
        F = f_penalty_bc(F,ndofn,nn,ns,sm,C,ndlc,dlcm)

    elif bc_type == 2:
        K = k_elim_bc()
        F = f_elim_bc()

    return K, F

def k_penalty_bc(Kmod,ndofn,ns,sm,C):

    """Apply BCs to Stiffness Matrix w/ Penalty Method

    Params:
        K (ndarray): global stiffness matrix w/o BCs
        ndof (int): number of DOFs per node for element type
        ns (int): number of supported nodes
        sm (ndarray): support matrix
        C (float): penalty scaling factor

    Returns:
        K (ndarray): global stiffness matrix with BCs applied
    """
    for i in range(ns):
        node = sm[i,0]
        sups = sm[i,1:]
        k = sups.shape[0]-1
        for n in range(sups.shape[0]):
            if sups[n] == 1:
                j = ndofn*node-k-1 #minus one converts to array index
                Kmod[j,j] = Kmod[j,j] + C
            k = k-1

    return Kmod

def f_penalty_bc(F,ndofn,nn,ns,sm,C,ndlc,dlcm):

    """Apply BCs to Load Vector w/ Penalty Method

    Params:
        F (ndarray): global load vector w/o BCs
        ndofn (int): number of DOFs per node for element type
        nn (int): number of nodes
        ns (int): number of supports
        sm (int): support matrix
        C (float): penalty scaling factor
        ndlc (int): number of support displacements
        dlcm (ndarray): support displacement matrix

    Returns:
        F (ndarray): global load vector w/ BCs

    Notes:
        1. Explanation of a: 
    """

    a = assemble_sup_disp_vector(ndofn,nn,ndlc,dlcm)

    for i in range(ns):
        node = sm[i,0]
        sups = sm[i,1:]
        k = sups.shape[0]-1
        for n in range(sups.shape[0]):
            if sups[n] == 1:
                j = ndofn*node-k-1 #minus one converts to array index
                F[j] = F[j] + C*a[j]
            k = k-1

    return F

def assemble_sup_disp_vector(ndofn,nn,ndlc,dlcm):

    """Creates support displacement vector

    Params:
        ndofn (int): number of DOFs per node for element type
        nn (int): number of nodes
        ndlc (int): number of support displacements
        dlcm (ndarray): support displacement matrix

    Returns:
        a (ndarray): support displacement vector

    """

    a = np.zeros((ndofn*nn))
    k = ndofn 

    for i in range(dlcm.shape[0]):
        dlm = dlcm[i,1:]
        node = int(dlm[0])
        dirc = int(dlm[1])
        mag = dlm[2]

        if dirc == 11:
            dof = ndofn*node - k
            a[dof] = mag
        if dirc == 12:
            dof = ndofn*node - (k-1)
            a[dof] = mag
        if dirc == 13:
            dof = ndofn*node - (k-2)
            a[dof] = mag
        if dirc == 21:
            dof = ndofn*node - (k-3)
            a[dof] = mag
        if dirc == 22:
            dof = ndofn*node - (k-4)
            a[dof] = mag
        if dirc == 23:
            dof = ndofn*node - (k-5)
            a[dof] = mag

    return a

def k_elim_bc():
    pass

def f_elim_bc():
    pass

