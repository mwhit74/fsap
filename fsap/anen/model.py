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
        """Generate structure coordinate vector."""
    
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
