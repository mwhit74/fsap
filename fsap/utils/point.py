"""
2D and 3D point classes. 
"""

#Need to consider the implications of refactoring the class to a 3D point
#Refactor class to 3D point only; set z as a named argument defaulting to 0

#Polar coordinate coversion?

class ThreeDPoint:
    """Class to represent a point in 3D space."""

    def __init__(self,x,y,z=0.0):
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

    def rect_2D(self):
        """Returns a tuple of x and y coordinates."""
        return (self.x, self.y)

    #def polar_2D(self):
    #    """Returns a tuple of magnitude and direction."""
    #    r = math.sqrt(math.pow(self.x,2) + math.pow(self.y,2))
    #    t = math.degrees(math.atan(self.x/self.y))

        return(r, t)

    def rect_3D(self):
        """Returns a tuple of x, y, and z coordinates."""
        return (self.x, self.y, self.z)

    #def polar_3D(self):
