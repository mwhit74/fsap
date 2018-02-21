from fsap.geom import line

class Element:
    def __init__(self, name, beg, end, matl, sp):
        self.name = name
        self.beg = beg
        self.end = end
        self.matl = matl
        self.sp = sp

        self.line = line.Line3D(self.beg, self.end)
        self.length = self.line.length()

    def tranformation_matrix():
        pass


    def local_stiffness_matrix():
        pass


    def global_stiffness_matrix():
        pass


    def local_member_end_forces():
        pass


    def global_member_end_forces():
        pass


    def local_member_end_displacements():
        pass


    def global_member_end_displacements():
        pass
