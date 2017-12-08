from fsap.utils.point import Point2D
import math
import functools
import pdb

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
        self.points = []
        self.convert_to_points(points)
        #self.start_pt = self.start_pt()
        self.centroid = self.centroid()
        self.order_points()
        self.add_end_pt()

        self.max_y = self.max_y()
        self.max_x = self.max_x()
        self.min_y = self.min_y()
        self.min_x = self.min_x()
        self.box = self.bounding_box()
        self.height = self.height()
        self.width = self.width()

        self.area = self.area()
        self.ena_x = self.ena_x()
        self.ena_y = self.ena_y()

        self.yt = self.yt()
        self.yb = self.yb()
        self.xl = self.xl()
        self.xr = self.xr()

        self.ixx_x = self.ixx_x()
        self.iyy_y = self.iyy_y()
        self.ixx_c = self.ixx_c()
        self.iyy_c = self.iyy_c()
        self.ixy_xy = self.ixy_xy()
        self.ixy_c = self.ixy_c()

        #self.sxxt = self.sxxt()
        #self.sxxb = self.sxxb()
        #self.syyl = self.syyl()
        #self.syyr = self.syyr()

        #self.rx = self.rx()
        #self.ry = self.ry()
        #self.rmin = self.rmin()


    def convert_to_points(self, points):
        """Converts list of tuples to list of Point2D objects"""
        if all(isinstance(pt, Point2D) for pt in points):
            self.points = points
        else:
            for pt in points:
                self.points.append(Point2D(pt[0], pt[1]))


    def start_pt(self):
        """Find point in lower left position"""
        l1 = sorted(self.points, key=lambda pt: pt.x)
        l2 = sorted(l1, key=lambda pt: pt.y)
        return l2[0]


    def order_points(self):
        """Order points in a counter-clockwise direction

        References:
        https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order
        https://en.wikipedia.org/wiki/Cross_product
        https://en.wikipedia.org/wiki/Shoelace_formula
        """
        self.points.sort(key=functools.cmp_to_key(self.less))


    def less(self,a, b):
        """Comparison function to determine the order of two points"""
        #if a is right of center and b is left of center, a is ccw from b
        if a.x - self.centroid.x >= 0.0 and b.x - self.centroid.x < 0.0:
            return 1
        #if a is left of center and b is right of center, a is cw from b
        if a.x - self.centroid.x < 0.0 and b.x - self.centroid.x >= 0:
            return -1
        #if a, b, and center lie on the same same vertical line
        if a.x - self.centroid.x == 0.0 and b.x - self.centroid.x == 0.0:
            #if a.y is greater than b.y, a.y is ccw from b.y otherwise a.y is cw
            #from b.y
            if a.y < b.y:
                return 1
            else:
                return -1

            #if not requiring points to be in first quadrant
            #if a.y - center.y >= 0 or b.y - center.y >=0:
            #    return a.y > b.y
            #return b.y > a.y

        #if a.x and b.x are on the same side of center calculate the
        #cross-product  of (center -> a) x (center -> b)
        det = ((a.x - self.centroid.x)*(b.y - self.centroid.y) - 
                (b.x - self.centroid.x)*(a.y - self.centroid.y))
        #if the cross-product is positive a is cw from b, otherwise a is ccw
        #from b
        if (det < 0):
            return 1
        else:
            return -1

        #a and b lie on the same line from the center
        #if a is farther than b, a is ccw from b, otherwise a is cw from b
        d1 = math.pow((a.x - self.centroid.x),2) + math.pow((a.y - self.centroid.y),2)
        d2 = math.pow((b.x - self.centroid.x),2) + math.pow((b.y - self.centroid.y),2)
        print d1 > d2
        if d1 > d2:
            return 1
        else:
            return -1


    def add_end_pt(self):
        self.points.append(self.points[0])


    def max_y(self):
        """Finds the maximum y-coordinate"""
        max_y_pt = max(self.points, key=lambda pt: pt.y)
        return max_y_pt.y


    def min_y(self):
        """Finds the minimum y-coordinate"""
        min_y_pt = min(self.points, key=lambda pt: pt.y)
        return min_y_pt.y


    def max_x(self):
        """Finds the maximum x-coordinate"""
        max_x_pt = max(self.points, key=lambda pt: pt.x)
        return max_x_pt.x


    def min_x(self):
        """Finds the minimum x-coordinate"""
        min_x_pt = min(self.points, key=lambda pt: pt.x)
        return min_x_pt.x


    def bounding_box(self):
        """Finds the max and min x and y coordinates"""
        return self.max_x, self.max_y, self.min_x, self.min_y


    def height(self):
        """Calculates the height of the polygon"""
        return self.max_y - self.min_y

    
    def width(self):
        """Calculates the width of the polygon"""
        return self.max_x - self.min_x


    def yt(self):
        """Calculates the distance from ENA to extreme top fibre"""
        return self.max_y - self.ena_y


    def yb(self):
        """Calculates the distance from ENA to extreme bottom fibre"""
        return self.ena_y - self.min_y


    def xr(self):
        """Calculates the distance from ENA to extreme right fibre"""
        return self.max_x - self.ena_x


    def xl(self):
        """Calculates the distance from ENA to extreme left fibre"""
        return self.ena_x - self.min_x


    def area(self):
        """Calculates the area of the polygon"""
        return -1*self.loop(self.area_eq)


    def area_eq(self, cur_x, cur_y, next_x, next_y):
        """Equation used to calculate area of polygon"""
        return (next_y - cur_y)*(next_x + cur_x)/2.0

            
    def ena_x(self):
        """Calculates elastic neutral axis in x-direction"""
        return -1.0/self.area*self.loop(self.ena_x_eq)
           

    def ena_x_eq(self, cur_x, cur_y, next_x, next_y):
        """Equation use to calculate elastic neutral axis in x-direction"""
        return ((next_y - cur_y)/8.0)*((next_x + cur_x)**2.0 + (next_x - cur_x)**2.0/3.0)
            

    def ena_y(self):
        """Calculates elastic neutral axis in y-direction"""
        return 1.0/self.area*self.loop(self.ena_y_eq)
        

    def ena_y_eq(self, cur_x, cur_y, next_x, next_y):
        return ((next_x - cur_x)/8.0)*((next_y + cur_y)**2.0 + (next_y - cur_y)**2.0/3.0)


    def centroid(self):
        """Calculates the centroid coordinate of the cross-section"""
        sum_x = 0.0
        sum_y = 0.0
        #subtract 1 from num points, last point closes polygon
        num_pts = len(self.points)
        for pt in self.points:
            sum_x = sum_x + pt.x
            sum_y = sum_y + pt.y

        xc = sum_x/num_pts
        yc = sum_y/num_pts

        return Point2D(xc, yc)
            
    def ixx_x(self):
        """Calculates second moment of area about x-axis (y=0)"""
        return self.loop(self.ixx_x_eq)
    

    def ixx_x_eq(self, cur_x, cur_y, next_x, next_y):
        """Equation used to calculate the second moment of area about x-axis"""
        return ((next_x - cur_x)*(next_y + cur_y)/24.0)*((next_y + cur_y)**2.0 + (next_y - cur_y)**2.0)
           

    def iyy_y(self):
        """Calculates second moment of area about y-axis (x=0)"""
        return -1*self.loop(self.iyy_y_eq)
           

    def iyy_y_eq(self, cur_x, cur_y, next_x, next_y):
        """Equation used to calculate the second moment of area about y-axis"""
        return ((next_y - cur_y)*(next_x + cur_x)/24.0)*((next_x + cur_x)**2.0 + (next_x - cur_x)**2.0)
            
            
    def ixx_c(self):
        """Calculates the second moment of area about x-axis of centroid"""
        return self.ixx_x - self.area*math.pow(self.ena_y,2)
            

    def iyy_c(self):
        """Calculates the second moment of area about y-axis of centroid"""
        return self.iyy_y - self.area*math.pow(self.ena_x,2)


    def ixy_xy(self):
        """Calculates the polar moment of area about the xy axis (0,0)"""
        return self.loop(self.ixy_xy_eq)


    def ixy_xy_eq(self, cur_x, cur_y, next_x, next_y):
        """Equation used to calc the polar moment area about the xy axis"""
        if next_x - cur_x == 0.0:
            a = 0.0
        else:
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
        """Calculates the polar moment of area about the centroid"""
        return self.ixy_xy + self.area*self.ena_x*self.ena_y


    def sxxt(self):
        """Calculates the elastic section modulus wrt extreme top fibre"""
        return self.ixx_c/self.yt


    def sxxb(self):
        """Calculates the elastic section modulus wrt extreme bottom fibre"""
        return self.ixx_c/self.yb


    def syyr(self):
        """Calculates the elastic section modulus wrt extreme right fibre"""
        return self.iyy_c/self.xr


    def syyl(self):
        """Calculates the elastic section modulus wrt extreme left fibre"""
        return self.iyy_c/self.xl


    def rx(self):
        """Calculates the radius of gyration about the x-axis of centroid"""
        return math.sqrt(self.ixx_c/self.area)


    def ry(self):
        """Calculates the radius of gyration about the y-axis of centroid"""
        return math.sqrt(self.iyy_c/self.area)


    def rmin(self):
        """Finds the minimum radius of gyration"""
        return min(self.rx, self.ry)


    def loop(self, func):
        """Loops thru points using given equation (func) to calc section prop"""
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
    points = [(1.0,1.0),(1.0,0.0),(0.0,1.0),(0.0, 0.0)]
    sp1 = SectionProperty(points)
    for pt in sp1.points:
        print pt

    print "\n"
    points = [(0.0, 0.0),(1.0,1.0),(1.0,0.0),(0.0,1.0)]
    sp2 = SectionProperty(points)
    for pt in sp2.points:
        print pt

    print "\n"
    points = [(1.0,0.0),(2.0,1.0),(3.0,0.0)]
    sp3 = SectionProperty(points)
    for pt in sp3.points:
        print pt

    print "\n"
    points = [(2.0,1.0),(3.0,0.0),(1.0,0.0)]
    sp4 = SectionProperty(points)
    for pt in sp4.points:
        print pt
