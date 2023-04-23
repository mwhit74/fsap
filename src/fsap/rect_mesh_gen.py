#domain width (x) and height (y)
lx = 100
ly = 10

#number of elements in x and y dirs
nx = 50
ny = 10

#element size
dx = lx/nx
dy = ly/ny

#aspect ratio
asr = max(dx,dy)/min(dx,dy) 

#number of steps in x and y
xs = lx/dx
yx = ly/dy

#number of global elements and nodes
ne = xs*ys
nn = (xs+1)*(ys+1)


def generate_nodes(nn,xs,ys)

    """Generates nodes for mesh

    Params:
        nn (int): total number of global nodes
        xs (int): number of divisions in the x-dir
        ys (int): number of divisions in the y-dir

    Returns:
        ncm (ndarray): nodal x and y coords

    """

    node = 1
    nc = np.empty((nn,3))
    for i in range(xs+1):
        for j in range(xs+1):
            xc = j*dx-dx #x-coord
            yc = i*dy-dy #y-coord
            ncm[n] = np.array([node,xc,yc])
            node = node + 1

    return ncm

def generate_element_connectivity(ne,ys,xs):

    """Generates nodes for mesh

    Params:
        ne (int): total number of global elements
        xs (int): number of divisions in the x-dir
        ys (int): number of divisions in the y-dir

    Returns:
        ecm (ndarray): element connectivity matrix

    """

    k = 1
    e = 1
    ecm = np.array((ne,5))
    for i in range(ys):
        for j in range(xs)
            n1 = k
            n2 = k+1
            n3 = k+(xs+2)
            n4 = k+(xs+1)
            ecm[e] = np.array([e,n1,n2,n3,n4])
            k = k + 1
            e = e + 1
        k = i*(xs+1)+1

    return ecm

def generate_element_property_matrix(ne,mpn,spn):

    """Generates element property assigment matrix

    Params:
        mpn (int): material property id number
        spn (int): section property id number

    Returns:
        epm (ndarray): element property assignment matrix
    """

    epm = np.array((ne,3))
    for e in range(ne):
        epm[e] = np.array([e,mpn,spn])

    return epm

def region_edge_nodes(nn,xs,ys):
    
    """Returns region edge nodes

    """
    be = np.empty(xs+1)
    te = np.empty(xs+1)
    re = np.empty(ys+1)
    le = np.empty(ys+1)

    nl = 1
    nh = nn
    for i in range(xs+1):
        be[i] = nl
        te[i] = nh
        nl+=1
        nh-=1

    ln = 1
    rn = xs+1
    for i in range(ys+1):
        le[i] = ln
        re[i] = rn
        ln = (i+1)*(xs+2)
        rn = (i+1)*(xs+1)

    return be, te, re, le
        
def region_corner_nodes(xs,ys):

    """Return region corner nodes

    """

    bl = 1
    br = xs+1
    tl = nn-xs-1
    tr = nn

    return bl,br,tl,tr

def plot_mesh():
    pass        
