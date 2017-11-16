from fsap.utils import Point2D
from math

#set up unit testing
#start working on docstrings for section property class

#additional functionality to add:
#   calc principle axes
#   calc PNA 
#       +create new class for composite x-section and use SectionProperty
#        as base class
#       +plastic moment of inertia(?)
#       +plastic section modulus
#   calc torsional properties
#       +polar moment of inertia
#       +open vs closed section
#   calc shear area
#       +will be different for different shapes
#       +may need to use first moment of inertia
#       +hollow tubes are different
#   first moment of area - Qx and Qy
#   shape factor

class SectionProperty:
    """
    Elastic section property class

    User Notes:
        1. All points must be located entirely within the first quadrant
        2. A group of points must be specified as a solid or void
    """
    def __init__(self, points):
        """Initialize section property class

        Args:
            points (list of tuples or list of Point2D): coordinate points
            representing the outline of the polygon
        Returns:
            None
        """
        if all(isinstance(points, Point2D):
            self.points = points
        else:
            self.points = convert_to_points(points)
       
        self.points = self.order_points()


    def convert_to_points(self):
        """Converts list of tuples to list of Point2D objects"""
        points = []
        for pt in points:
            points.append(Point2D(pt[0], pt[1]))
        return points


    def order_points(self):
        https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order
        https://en.wikipedia.org/wiki/Cross_product


    def max_y(self):
        max_y = None
        for pt in self.points:
            if pt.y > max_y:

        return max_y


    def min_y(self):
        min_y = None
        for pt in self.points:
            if pt.y < min_y:
                min_y = pt.y


    def max_x(self):
        max_x = None
        for pt in self.points:
            if pt.x > max_x:
                max_x = pt.x


    def min_x(self):
        min_x = None
        for pt in self.points:
            if pt.x < min_x:
                min_x = pt.x


    def bounding_box(self):
        """Finds the max and min x and y coordinates"""
        return self.max_x(), self.max_y(), self.min_x(), self.min_y()


    def height(self):
        """Calculates the height of the polygon"""
        return self.max_y() - self.min_y()

    
    def width(self):
        """Calculates the width of the polygon"""
        return self.max_x() - self.min_x()


    def yt(self):
        """Calculates the distance from ENA to extreme top fibre"""
        return self.max_y() - self.ena_y()


    def yb(self):
        """Calculates the distance from ENA to extreme bottom fibre"""
        return self.eny_y() - self.min_y()


    def xr(self):
        """Calculates the distance from ENA to extreme right fibre"""
        return self.max_x() - self.ena_x()


    def xl(self):
        """Calculates the distance from ENA to extreme left fibre"""
        return self.ena_x() - self.min_x()


    def area(self):
        """Calculates the area of the polygon"""
        return -1*self.loop(self.area_eq())


    def area_eq(self, cur_x, cur_y, next_x, next_y):
        """Equation used to calculate area of polygon"""
        return (next_y - cur_y)*(next_x + cur_x)/2.0

            
    def ena_x(self):
        """Calculates elastic neutral axis in x-direction"""
        return -1.0/self.area()*self.loop(self.ena_x_eq())
           

    def ena_x_eq(self, cur_x, cur_ye, next_x, next_y):
        """Equation use to calculate elastic neutral axis in x-direction"""
        return ((next_y - cur_y)/8.0)*((next_x + cur_x)**2.0 + (next_x - cur_x)**2.0/3.0)
            

    def ena_y(self):
        """Calculates elastic neutral axis in y-direction"""
        return 1.0/self.area()*self.loop(self.ena_y_eq())
        

    def ena_y_eq(self, cur_x, cur_y, next_x, next_y):
        return ((next_x - cur_x)/8.0)*((next_y + cur_y)**2.0 + (next_y - cur_y)**2.0/3.0)


    def centroid(self):
        return Point2D(self.ena_x(), self.ena_y())
            
    def ixx_x(self):
        """Calculates second moment of area about x-axis (y=0)"""
        return self.loop(self.ixx_x_eq())
    

    def ixx_x_eq(self, cur_x, cur_y, next_x, next_y):
        return ((next_x - cur_x)*(next_y + cur_y)/24.0)*((next_y + cur_y)**2.0 + (next_y - cur_y)**2.0)
           

    def iyy_y(self):
        return -1*self.loop(self.iyy_y_eq())
           

    def iyy_y_eq(self, cur_x, cur_y, next_x, next_y):
        return ((next_y - cur_y)*(next_x + cur_x)/24.0)*((next_x + cur_x)**2.0 + (next_x - cur_x)**2.0)
            
            
    def ixx_c(self):
        return self.ixx_x() - self.area()*math.pow(self.ena_y(),2)
            

    def iyy_c(self):
        return self.iyy_y() - self.area()*math.pow(self.ena_x(),2)


    def ixy_xy(self):
        return self.loop(self.ixy_xy_eq())


    def ixy_xy_eq(self):
        a = 1.0/(next_x - cur_x)
        b = 1.0/8.0
        c = math.pow(next_y - cur_y,2)
        d = next_x + cur_x
        e = math.pow(next_x,2) + math.pow(cur_x,2)
        f = 1.0/3.0
        g = next_y - cur_y
        h = next_x*cur_y - cur_x*next_y
        i = math.pow(next_x,2) + next_x*cur_x + math.pow(cur_x,2)
        j = 1.0/4.0
        k = math.pow(next_x*cur_y - cur_x*next_y,2)
        l = next_x + cur_x

        return a*((b*c*d*e)+(f*g*h*i)+(j*k*l))

    def ixy_c(self):
        return self.ixy_xy() + self.area()*self.ena_x()*self.ena_y()


    def sxxt(self):
        return self.ixx_c()/self.yt()


    def sxxb(self):
        return self.ixx_c()/self.yb()


    def syyr(self):
        return self.iyy_c()/self.xr()


    def syyl(self):
        return self.iyy_c()/self.xl()


    def rx(self):
        return math.sqrt(self.ixx_c()/self.area())


    def ry(self):
        return math.sqrt(self.iyy_c()/self.area())


    def rmin(self):
        return min(self.rx(), self.ry())


    def loop(self, func):
            var = 0.0
           
            cur_pt = self.points[0]
            cur_x = cur_pt.x
            cur_y = cur_pt.y
            
            for pt in self.points:
                    next_x = pt.x
                    next_y = pt.y
                    
                    var = var + func(cur_x, cur_y, next_x, next_y)
                            
                    cur_x = next_x
                    cur_y = next_y
            
            return var

	
if __name__ == "__main__":
	#xy_coords = [[0.0,0.0],[2.0,0.0],[2.0,2.0],[0.0,2.0],[0.0,0.0]] #rectangle
	#xy_coords = [[0,0],[2,0],[1,2],[0,0]] #triangle
	xy_coords = [[0.0,0.0],[7.0,0.0],[8.26,3.5],[0.0, 3.5],[0.0,0.0]]
	
	print "Area: " + str(area(xy_coords))
	print "ena in x: " + str(ena_x(xy_coords))
	print "ena in y: " + str(ena_y(xy_coords))
	print "I about x: " + str(I_about_x(xy_coords))
	print "I about y: " + str(I_about_y(xy_coords))
	print "I about x at centroid: " + str(I_about_x_centroid(xy_coords))
	print "I about y at centroid: " + str(I_about_y_centroid(xy_coords))
