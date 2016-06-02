class XYPoint:
    def __init__(self,x,y):
        self.x = x
        self.y = y

    def coord(self):
        return (self.x,self.y)

class XYZPoint:
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z

    def coord(self):
        return (self.x, self.y, self.z)


