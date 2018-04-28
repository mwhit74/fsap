class Joint(object):
    def __init__(name, x, y, sup, z=0.0)
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.sup = sup
        self.scv = None

class TrussJoint(Joint):
    def __init__(self, x, y, dx, dy):
        sup = [dx, dy]
        Joint.__init__(x, y, sup)

    def set_scv(self, sc_dx, sc_dy):
        self.scv = (sc_dx, sc_dy)

class TwoDFrameJoint(Joint):
    def __init__(self, x, y, dx, dy, rz):
        sup = [dx, dy, rz]
        Joint.__init__(x, y, sup)

    def set_scv(self, sc_dx, sc_dy, sc_rz):
        self.scv = (sc_dx, sc_dy, sc_rz)

class ThreeDFrameJoint(Joint):
    def __init__(self, x, y, z, dx, dy, dz, rx, ry, rz):
        sup = [dx, dy, dz, rx, ry, rz]
        Joint.__init__(x, y, sup, z=z)
        
    def set_scv(self, sc_dx, sc_dy, sc_dz, sc_rx, sc_ry, sc_rz):
        self.scv = (sc_dx, sc_dy, sc_dz, sc_rx, sc_ry, sc_rz)
