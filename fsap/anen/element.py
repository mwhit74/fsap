import math

class Element:
    def __init__(self, name, beg_jt, end_jt, sepr, matl):
        self.name = name
        self.beg_jt = beg_jt
        self.end_jt = end_jt
        self.sepr = sepr
        self.matl = matl
        self.k = None

        self.length = self.length()

        self.alpha = self.alpha()
        self.beta = self.beta()
        self.gamma = self.gamma

    def length(self):
        return math.sqrt(math.pow(self.beg_jt.x - self.end_jt.x,2) +
                         math.pow(self.beg_jt.y - self.end_jt.y,2) +
                         math.pow(self.beg_jt.z - self.end_jt.z,2))

    def alpah(self):
        pass

    def beta(self):
        pass

    def gamma(self):
        pass

    def scv(self):
        return join(self.beg_jt.scv,self.end_jt.scv)

class TrussElement(Element):
    def __init__(self, name, beg_jt, end_jt, sepr, matl):
        Element.__init__(name, beg_jt, end_jt, sepr, matl)

    def assemble_transformation_matrix(self):
        pass

    def assemble_local_stiffness_matrix(self):
        pass

class TwoDFrameElement(Element):
    def __init__(self, name, beg_jt, end_jt, sepr, matl):
        Element.__init__(name, beg_jt, end_jt, sepr, matl)

    def assemble_transformation_matrix(self):
        pass

    def assemble_local_stiffness_matrix(self):
        pass

class ThreeDFrameElement(Element):
    def __init__(self, name, beg_jt, end_jt, sepr, matl):
        Element.__init__(name, beg_jt, end_jt, sepr, matl)

    def assemble_transformation_matrix(self):
        pass

    def assemble_local_stiffness_matrix(self):
        pass
