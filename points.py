"""
XY and XYZ point classes. 
"""

__author__ = "MLW"
__copyright__ = ""
__credits__ = ["MLW"]

__license__ = ""
__version__ = "1.0"
__maintainer__ = "MLW"
__email__ = "mwhit74@gmail.com"
__status__ = "Development"


class XYPoint:
    """Class to represent a point in 2D space."""

    def __init__(self,x,y):
        """Initialize the class arguments x and y

        Args:
            x (float): x-coordinate
            y (float): y-coordinate
        """
        
        self.x = x
        self.y = y

    def coord(self):
        """Returns a tuple of x and y coordinates."""
        return (self.x,self.y)

class XYZPoint:
    """Class to represent a point in 3D space."

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

    def coord(self):
        """Returns a tuple of x, y, and z coordinates."""
        return (self.x, self.y, self.z)


