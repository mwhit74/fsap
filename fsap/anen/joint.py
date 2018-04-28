import fsap.point as pt

class Joint:
    def __init__(name, x, y, z=0.0, ux, uy, uz=1, rx=1, ry=1, rz=1):
        self.name = name
        self.coord = pt.Point3D(x,y,z)
        self.ux = ux
        self.uy = uy
        self.uz = uz
        self.rx = rx
        self.ry = ry
        self.rz = rz

