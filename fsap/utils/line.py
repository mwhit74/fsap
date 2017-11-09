import math

"
Line class
"

class Line2D:
    "Class to represent a line in space."

    def __init__(self, pt1, pt2):
        "Initialize the class arguments pt1 and pt2.

        Args:
            pt1 (Point2D): first point
            pt2 (Point2D): second point

        *See point.py for Point2D class definition. 
        """

        self.pt1 = pt1
        self.pt2 = pt2

    def __str__(self):
        """String representation of a line."""
        return "".format(self.pt1.x, self.pt1.y, self.pt2.x, self.pt2.y)

    def length(self):
        """Calculate the length of a line."""
        return math.sqrt(pow(self.pt1.x - self.pt2.x, 2) +\
                pow(self.pt1.y-self.pt2.y,2))


class Line3D:
    "Class to represent a line in space."

    def __init__(self, pt1, pt2):
        "Initialize the class arguments pt1 and pt2.

        Args:
            pt1 (Point3D): first point
            pt2 (Point3D): second point

        *See point.py for Point3D class definition. 
        """

        self.pt1 = pt1
        self.pt2 = pt2

    def __str__(self):
        """String representation of a line."""
        return "".format(self.pt1.x, self.pt1.y, self.pt1.z,
                         self.pt2.x, self.pt2.y, self.pt2.z)

    def length(self):
        """Calculate the length of a line."""
        return math.sqrt(pow(self.pt1.x - self.pt2.x,2) +
                         pow(self.pt1.y - self.pt2.y,2) +
                         pow(self.pt1.z - self.pt2.z,2))



