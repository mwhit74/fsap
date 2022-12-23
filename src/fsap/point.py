import numpy as np

"""
3D point class. 
"""

#Polar coordinate coversion?

class Point3D:
    """Class to represent a point in 3D space."""

    def __init__(self,x,y,z = 0.0):
        """Initialize the class arguments x, y, and z.

        Args:
            x (float): x-coordinate
            y (float): y-coordinate
            z (float): z-coordinate; defaults to 0.0
        """

        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        """String representation of a 3D point."""
        return "({0:.3f},{1:.3f},{2:.3f})".format(self.x, self.y, self.z)

    def cc(self):
        """Returns a tuple of x, y, and z carteisan coordinates."""
        return np.array([self.x, self.y, self.z])


