import math
import numpy as np

class Element:
    def __init__(self, name, beg_jt, end_jt, sepr, matl):
        self.name = name
        self.beg_jt = beg_jt
        self.end_jt = end_jt
        self.sepr = sepr
        self.matl = matl

        self.length = self.calc_length()
        self.uv = self.calc_unit_vector()
        self.scv = self.calc_scv()
        
        self.t = None
        self.k = None

    def calc_length(self):
        return math.sqrt(math.pow(self.beg_jt.x - self.end_jt.x,2) +
                         math.pow(self.beg_jt.y - self.end_jt.y,2) +
                         math.pow(self.beg_jt.z - self.end_jt.z,2))
        
    def calc_unit_vector(self):
        return np.array([(self.end_jt.x - self.beg_jt.x)/self.length,
                         (self.end_jt.y - self.end_jt.y)/self.length,
                         (self.end_jt.z - self.end_jt.z)/self.length])

    def calc_scv(self):
        return join(self.beg_jt.scv,self.end_jt.scv)

class BarElement(Element):
    def __init__(self, name, beg_jt, end_jt, sepr, matl):
        Element.__init__(name, beg_jt, end_jt, sepr, matl)

    def assemble_transformation_matrix(self):
        self.t = np.array([[self.uv[0], self.uv[1], self.uv[2], 0., 0., 0.],
                           [0., 0., 0., self.uv[0], self.uv[1], self.uv[2]]])
        
    def assemble_local_stiffness_matrix(self):
        p = self.matl.e*self.sepr.area/self.length
        
        k = np.array([1., -1.],[-1., 1.])
        
        self.k = np.multiply(p, k)
    
class ThreeDThickBeamElement(Element):
    def __init__(self, name, beg_jt, end_jt, sepr, matl):
        Element.__init__(name, beg_jt, end_jt, sepr, matl)
        
        self.alpha = self.calc_alpha()
        self.beta = self.calc_beta()
        self.gamma = self.calc_gamma()
        
    def alpah(self):
        pass

    def beta(self):
        pass

    def gamma(self):
        pass

    def assemble_transformation_matrix(self):
        pass

    def assemble_local_stiffness_matrix(self):
        pass

class ThreeDThinBeamElement(Element):
    def __init__(self, name, beg_jt, end_jt, sepr, matl):
        Element.__init__(name, beg_jt, end_jt, sepr, matl)

    def assemble_transformation_matrix(self):
        pass

    def assemble_local_stiffness_matrix(self):
        pass
