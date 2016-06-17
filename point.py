"""
XY and XYZ point classes. 
"""

class XYPoint:
    """Class to represent a point in 2D space."""

    def __init__(self,x,y):
        """Initialize the class arguments x and y.

        Args:
            x (float): x-coordinate
            y (float): y-coordinate
        """
        
        self.x = x
        self.y = y

    def __str__(self):
        """String representation of a 2D point."""
        return "({0:.3f}, {1:.3f})".format(self.x,self.y)

    def coord(self):
        """Returns a tuple of x and y coordinates."""
        return (self.x,self.y)

class XYZPoint:
    """Class to represent a point in 3D space."""

    def __init__(self,x,y,z):
        """Initialize the class arguments x, y, and z.

        Args:
            x (float): x-coordinate
            y (float): y-coordinate
            z (float): z-coordinate
        """

        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        """String representation of a 3D point."""
        return "({0:.3f}, {1:.3f}, {2:.3f})".format(self.x, self.y, self.z)

    def coord(self):
        """Returns a tuple of x, y, and z coordinates."""
        return (self.x, self.y, self.z)


