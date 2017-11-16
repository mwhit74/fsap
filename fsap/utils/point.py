import numpy as np

"""
2D and 3D point classes. 
"""

#Need to consider the implications of refactoring the class to a 3D point
#Refactor class to 3D point only; set z as a named argument defaulting to 0

#Polar coordinate coversion?

class Point2D:
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
        return "{0:.3f} {1:.3f}".format(self.x,self.y)

    def cartesian_coord(self):
        """Returns a tuple of x and y cartesian coordinates."""
        return np.array([self.x,self.y])
    
class Point3D:
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
        return "{0:.3f} {1:.3f} {2:.3f}".format(self.x, self.y, self.z)

    def coord(self):
        """Returns a tuple of x, y, and z coordinates."""
        return np.array([self.x, self.y, self.z])


