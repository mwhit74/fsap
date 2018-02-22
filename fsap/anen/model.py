import numpy as np

class Model:
    def __init__(self, jts, sup, matl, sect, elem, load):
        self.jts = jts
        self.sup = sup
        self.matl = matl
        self.sect = sect
        self.elem = elem
        self.load = load

    def run():
        pass
        
    def num_dof(self, sup, num_scj, num_jt):
        """Calculate the number of degrees of freedom of the structure.
        
        Args:
            sup (numpy array): 2D array of joint number, x support condition, y
                              support condition
            num_scj (int): number of structure coordinated per joint; number of
                           degrees of freedom per joint
            num_jt (int): number of joints in the structure
    
        Returns:
            ndof (int): number of degrees of freedom of the structure
    
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
    
    def scv(self, sup, num_scj, num_jt, num_dof):
        """Generate structure coordinate vector.
        
        Args:
            sup (list): 2D array of joint number, x support condition, y support
                        condition
            num_scj (int): number of structure coordinates per joint; number of
                           degrees of freedom per joint
            num_jt (int): number of joints in the structure
            num_dof (int): number of degrees of freedom of the structure
    
        Returns:
            str_cv (numpy array): 1D array (vector) of the numbered DOFs for
                                  the entire structure to the joints
        """
    
        str_cv = np.empty([num_jt*num_sj],dtype=int)
        j = 0
        k = num_dof
        for i in range(num_jt):
            if i in sup[0]:
                for x in range(num_scj):
                    y = (i - 1)*num_scj + x #location of current dof in scv
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

    
    def assemble_stiffness(num_dof, num_scj):
        """Assemble structure stiffness matrix."""
    
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

