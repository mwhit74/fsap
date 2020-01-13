from fsap.geom.point import Point3D
from fsap.anen.suppoprt import Support

class Joint(Point3D):
    def __init__(self, name, x, y, z, sup):
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.sup = sup
        self.scv = None



