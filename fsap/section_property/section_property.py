from fsap.utils import Point2D
from math import sqrt


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


    """
    def __init__(self, points):
        if all(isinstance(points, Point2D):
            self.points = points
        else:
            self.points = convert_to_points(points)
       
        self.order_points = self.order_points()

        self.max_y, self.min_y, self.max_x, self.min_x = self.bounding_box()

        self.height = self.height()
        self.width = self.width()

        self.area = self.area()
        self.ena_x = self.ena_x()
        self.ena_y = self.ena_y()
        self.ixx_x = self.i_about_x()
        self.iyy_y = self.i_about_y()
        self.ixx_c = self.i_about_x_centroid()
        self.iyy_c = self.i_about_y_centroid()
        self.ixy_xy = self.ixy_about_xy()
        self.ixy_c = self.ixy_about_centroid()

        self.yt = self.yt()
        self.yb = self.yb()
        self.xr = self.xr()
        self.xl = self.xl()

        self.sxxt = self.ixx_c/self.yt
        self.sxxb = self.ixx_c/self.yb
        self.syyr = self.iyy_c/self.xr
        self.syyl = self.iyy_c/self.xl

        self.rx = math.sqrt(self.ixx_c/self.area)
        self.ry = math.sqrt(self.iyy_c/self.area)
        self.rmin = min(self.rx, self.ry)

    def convert_to_points(self):
        points = []
        for pt in points:
            points.append(Point2D(pt[0], pt[1]))
        return points


    def bounding_box(self):
        max_y = 0.0
        min_y = 0.0
        max_x = 0.0
        min_x = 0.0
        for pt in self.points:
            if pt.y > max_y:
                max_y = pt.y
            if pt.y < min_y:
                min_y = pt.y
            if pt.x > max_x:
                max_x = pt.x
            if pt.y < min_y:
                min_x = pt.x

        return max_y, min_y, max_x, min_x


    def height(self):
        return self.max_y - self.min_y

    
    def width(self):
        return self.max_x - self.min_x

    #not sure if distance from centroid to extreme points will work correctly
    #with respect the signs working right
    #do these need to be with respect to the centoid of the section to the get
    #signs on the section modulus correct?
    def yt(self):
        return self.max_y - self.ena_y


    def yb(self):
        return self.eny_y - self.min_y


    def xr(self):
        return self.max_x - self.ena_x


    def xl(self):
        return self.ena_x - self.min_x


    def area(self):
        return -1*self.loop(self.area_eq)
           

    def area_eq(self, cur_x, cur_y, next_x, next_y):
        return (next_y - cur_y)*(next_x + cur_x)/2.0

            
    def ena_x(self):
        return -1.0/self.area*self.loop(self.ena_x_eq)
           

    def ena_x_eq(self, cur_x, cur_y, next_x, next_y):
        return ((next_y - cur_y)/8.0)*((next_x + cur_x)**2.0 + (next_x - cur_x)**2.0/3.0)
            

    def ena_y(self):
        return 1.0/self.area*self.loop(self.ena_y_eq)
        

    def ena_y_eq(self, cur_x, cur_y, next_x, next_y):
        return ((next_x - cur_x)/8.0)*((next_y + cur_y)**2.0 + (next_y - cur_y)**2.0/3.0)
            
            
    def i_about_x(self):
        return self.loop(self.i_about_x_eq)
    

    def i_about_x_eq(self, cur_x, cur_y, next_x, next_y):
        return ((next_x - cur_x)*(next_y + cur_y)/24.0)*((next_y + cur_y)**2.0 + (next_y - cur_y)**2.0)
           

    def i_about_y(self):
        return -1*self.loop(self.i_about_y_eq)
           

    def i_about_y_eq(self, cur_x, cur_y, next_x, next_y):
        return ((next_y - cur_y)*(next_x + cur_x)/24.0)*((next_x + cur_x)**2.0 + (next_x - cur_x)**2.0)
            
            
    def i_about_x_centroid(self):
        return self.i_about_x - self.area*math.pow(self.ena_y,2)
            

    def i_about_y_centroid(self):
        return self.i_about_y - self.area*math.pow(self.ena_x,2)


    def ixy_about_xy(self):
        return self.loop(self.ixy_about_xy_eq)


    def ixy_about_xy_eq(self):
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

    def ixy_about_centroid(self):
        return self.ixy_about_xy + self.area*self.ena_x*self.ena_y

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

    #def graham_scan(xy_coords):
    #        x = xy_coords[0][0]
    #        y = xy_coords[0][1]

    #        for xy in xy_coords:
    #                if xy[1] < y:
    #                        y = xy[1]
    #                        x = xy[0]
    #                if xy[1] = y and xy[0] < x:
    #                        y = xy[1]
    #                        x = xy[0]

    #def ccw(p1, p2, p3):

	
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
