from fsap.geom.point import Point3D
import math
import functools
import pdb
import matplotlib.pyplot as plt

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
#   voids
#   circles

class SP:
    """
    Elastic section property class

    User Notes:
        1. All points must be located entirely within the first quadrant
    """
    def __init__(self, points):
        """Initialize section property class
        
        Args:
            points (list of tuples or list of Point3D): coordinate points
                representing the outline of the polygon
        Returns:
            None
        Notes:
        1. Points for a solid polygon must be entered in a clockwise direction.
        2. Points entereed for a void space must be entered in a
           counter-clockwise direction.
        3. The points entered as a list of tuples can be entered as a
           point on a 2D coordinate system, e.g. (x,y). The z-coordinate
           can be truncated and the Point3D class will automatically add
           a z-coordinate of 0.0.
        """
        self.points = []
        self.convert_to_points(points)
        self.centroid = self.centroid()
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

        self.sxxt = self.sxxt()
        self.sxxb = self.sxxb()
        self.syyl = self.syyl()
        self.syyr = self.syyr()

        self.rx = self.rx()
        self.ry = self.ry()
        self.rmin = self.rmin()


    def convert_to_points(self, points):
        """Converts list of tuples to list of Point3D objects"""
        if all(isinstance(pt, Point3D) for pt in points):
            self.points = points
        else:
            for pt in points:
                self.points.append(Point3D(pt[0], pt[1]))


    def order_points(self):
        """Order points in a clockwise direction for a convex polygon

        References:
        https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order
        https://en.wikipedia.org/wiki/Cross_product
        https://en.wikipedia.org/wiki/Shoelace_formula
        """
        self.points.sort(key=functools.cmp_to_key(self.less))


    def less(self,a, b):
        """Comparison function to determine the order of two points"""
        #if a is left of center and b is right of center, a is cw from b
        if b.x - self.centroid.x >= 0.0 and a.x - self.centroid.x < 0.0:
            return 1
        #if a is right of center and b is left of center, a is ccw from b
        if b.x - self.centroid.x < 0.0 and a.x - self.centroid.x >= 0:
            return -1
        #if a, b, and center lie on the same same vertical line
        if b.x - self.centroid.x == 0.0 and a.x - self.centroid.x == 0.0:
            #if a.y is greater than b.y, a.y is cw from b.y otherwise a.y is ccw
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
        if (det > 0):
            return 1
        else:
            return -1

        #a and b lie on the same line from the center
        #if a is farther than b, a is ccw from b, otherwise a is cw from b
        d1 = math.pow((a.x - self.centroid.x),2) + math.pow((a.y - self.centroid.y),2)
        d2 = math.pow((b.x - self.centroid.x),2) + math.pow((b.y - self.centroid.y),2)
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

        return Point3D(xc, yc)
            
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

    def all_props(self):
        """Returns a dictonary of all properties."""
        props = {}

        props["points"] = [pt.cc() for pt in self.points]

        props["max_y"] = self.max_y
        props["max_x"] = self.max_x
        props["min_y"] = self.min_y
        props["min_x"] = self.min_x
        props["bounding_box"] = self.box
        props["height"] = self.height
        props["width"] = self.width

        props["area"] = self.area
        props["ena_x"] = self.ena_x
        props["ena_y"] = self.ena_y

        props["yt"] = self.yt
        props["yb"] = self.yb
        props["xl"] = self.xl
        props["xr"] = self.xr

        props["ixx_x"] = self.ixx_x
        props["iyy_y"] = self.iyy_y
        props["ixx_c"] = self.ixx_c
        props["iyy_c"] = self.iyy_c
        props["ixy_xy"] = self.ixy_xy
        props["ixy_c"] = self.ixy_c

        props["sxxt"] = self.sxxt
        props["sxxb"] = self.sxxb
        props["syyl"] = self.syyl
        props["syyr"] = self.syyr

        props["rx"] = self.rx
        props["ry"] = self.ry
        props["rmin"] = self.rmin

        return props

    def plot_section(self, size_x = 10, size_y = 10):
        """Plots the section graphically."""
        
        #matplotlib patches and paths for plotting enclosed shapes
        
        pts = [pt.cc() for pt in self.points]
        xc = [pt[0] for pt in pts]
        yc = [pt[1] for pt in pts]
        
        plt.plot(xc, yc, "o")

        plt.xlabel("Width")
        plt.ylabel("Height")
        plt.title("Section Plot")
        plt.legend()
        plt.rcParams["figure.figsize"] = (size_x, size_y)
        

    def __str__(self):
        """Provides summary of all properties."""
        props = self.all_props()

        out_str = ""

        for p in props.keys():
            out_str = out_str + p + ": " + str(props[p]) + "\n"

        return out_str
    
class ISP(SP):
    
    def __init__(self, bf_top, tf_top, t_web, d_web, 
                 bf_bott=None, tf_bott=None,
                 bf_top_cp=0.0, tf_top_cp=0.0, 
                 bf_bott_cp=0.0, tf_bott_cp=0.0):
        
        self.bf_top = bf_top
        self.tf_top = tf_top
        self.t_web = t_web
        self.d_web = d_web
        self.bf_bott = bf_bott
        self.tf_bott = tf_bott
        
        self.bf_top_cp = bf_top_cp
        self.tf_top_cp = tf_top_cp
        self.bf_bott_cp = bf_bott_cp
        self.tf_bott_cp = tf_bott_cp
        
        self.pts = []
        self.calc_points()
        
        SP.__init__(self, self.pts)
        
    @property
    def bf_bott(self):
        return self.__bf_bott
    
    @bf_bott.setter
    def bf_bott(self, bf_bott):
        if bf_bott == None:
            self.__bf_bott = self.bf_top
        else:
            self.__bf_bott = bf_bott

    @property
    def tf_bott(self):
        return self.__tf_bott
    
    @tf_bott.setter
    def tf_bott(self, tf_bott):
        if tf_bott == None:
            self.__tf_bott = self.tf_top
        else:
            self.__tf_bott = tf_bott            
    
        
    def calc_points(self):
        
        if self.bf_top_cp >= self.bf_bott_cp:  
            
            d1 = self.bf_top_cp/2. - self.bf_bott_cp/2.
            d2 = self.bf_top_cp/2. - self.bf_bott/2.
            d3 = self.tf_bott_cp + self.tf_bott
            d4 = self.bf_top_cp/2. - self.t_web/2.
            d5 = d3 + self.d_web
            d6 = self.bf_top_cp/2. - self.bf_top/2.
            d7 = d5 + self.tf_top
            d8 = d7 + self.tf_top_cp
            d9 = self.bf_top_cp/2. + self.bf_top/2.
            d10 = self.bf_top_cp/2. + self.t_web/2.
            d11 = self.bf_top_cp/2. + self.bf_top/2.
            d12 = self.bf_top_cp/2. + self.bf_bott_cp/2.

            
            self.pts.append((d1, 0.))
            self.pts.append((d1, self.tf_bott_cp))
            self.pts.append((d2, self.tf_bott_cp))
            self.pts.append((d2, d3))
            self.pts.append((d4, d3))
            self.pts.append((d4, d5))
            self.pts.append((d6, d5))
            self.pts.append((d6, d7))
            self.pts.append((0., d7))
            self.pts.append((0., d8))
            self.pts.append((self.bf_top_cp, d8))
            self.pts.append((self.bf_top_cp, d7))
            self.pts.append((d9, d7))
            self.pts.append((d9, d5))
            self.pts.append((d10, d5))
            self.pts.append((d10, d3))
            self.pts.append((d11, d3))
            self.pts.append((d11, self.tf_bott_cp))
            self.pts.append((d12, self.tf_bott_cp))
            self.pts.append((d12, 0.))
            
            
        elif self.bf_top_cp < self.bf_bott_cp:
            print("Not defined")
        
    def __str__(self):
        """Provides summary of all properties."""
        
        user_input = (f'bf_top_cp = {self.bf_top_cp:.2f}, ' +
                      f'tf_top_cp = {self.tf_top_cp:.2f}, ' +
                      f'bf_top = {self.bf_top:.2f}, ' +
                      f'tf_top = {self.tf_top:.2f}, ' +
                      f't_web = {self.t_web:.2f}, ' +
                      f'd_web = {self.d_web:.2f}, ' +
                      f'bf_bott = {self.bf_bott:.2f}, ' +
                      f'tf_bott = {self.tf_bott:.2f}, ' +
                      f'bf_bott_cp = {self.bf_bott_cp:.2f}, ' +
                      f'tf_bott_cp = {self.tf_bott_cp:.2f} \n\n')
        
        all_props = super().__str__()
        
        return user_input + all_props
    
class TSP(SP):
    
    def __init__(self, bf, tf, t_stem, d_stem, bf_cp = 0.0, tf_cp = 0.0):
        
        self.bf_cp = bf_cp
        self.tf_cp = tf_cp
        self.bf = bf
        self.tf = tf
        self.t_stem = t_stem
        self.d_stem = d_stem

        
        self.pts = []
        self.calc_points()
        
        SP.__init__(self, self.pts)
        
    def calc_points(self):
        
        d1 = self.bf_cp/2. - self.tf_cp/2.
        d2 = self.bf_cp/2. - self.bf/2.
        d3 = self.d_stem + self.tf
        d4 = d3 + self.tf_cp
        d5 = self.bf_cp/2. + self.bf/2.
        d6 = self.bf_cp/2. + self.tf_cp/2.
        
        self.pts.append((d1, 0.))
        self.pts.append((d1, self.d_stem))
        self.pts.append((d2, self.d_stem))
        self.pts.append((d2, d3))
        self.pts.append((0., d3))
        self.pts.append((0., d4))
        self.pts.append((self.bf_cp, d4))
        self.pts.append((self.bf_cp, d3))
        self.pts.append((d5, d3))
        self.pts.append((d5, self.d_stem))
        self.pts.append((d6, self.d_stem))
        self.pts.append((d6, 0.))
        
    def __str__(self):
        """Provides summary of all properties."""
        
        user_input = (f'bf_cp = {self.bf_cp:.2f},' +
                      f'tf_cp = {self.tf_cp:.2f},' +
                      f'bf = {self.bf:.2f}, ' +
                      f'tf = {self.tf:.2f}, ' +
                      f't_stem = {self.t_stem:.2f}, ' +
                      f'd_stem = {self.d_stem:.2f} \n\n')
        
        all_props = super().__str__()
        
        return user_input + all_props
    
class RectSP(SP):
    def __init__(self, b, t, b_cp = 0.0, t_cp = 0.0):
        
        self.b_cp = b_cp
        self.t_cp = t_cp
        self.b = b
        self.t = t
        
        self.pts = []
        self.calc_points()
        
        SP.__init__(self, self.pts)
        
    def calc_points(self):
        
        d1 = self.b_cp/2. - self.b/2.
        d2 = self.t + self.t_cp
        d3 = self.b_cp/2. + self.b/2.
        
        self.pts.append((d1, 0.))
        self.pts.append((d1, self.t))
        self.pts.append((0., self.t))
        self.pts.append((0., d2))
        self.pts.append((self.b_cp, d2))
        self.pts.append((self.b_cp, self.t))
        self.pts.append((d3, self.t))
        self.pts.append((d3, 0.))
        
    def __str__(self):
        """Provides summary of all properties."""
        
        user_input = (f'b_cp = {self.b_cp:.2f},' +
                      f't_cp = {self.t_cp:.2f},' +
                      f'b = {self.b:.2f}, ' +
                      f't = {self.t:.2f} \n\n')
        
        all_props = super().__str__()
        
        return user_input + all_props
    
class Box(SP):
    
    def __init__(self, bf_top, tf_top, t_web, d_web, bf_bott = None, 
                 tf_bott = None, bf_top_cp = 0.0, tf_top_cp = 0.0,
                 bf_bott_cp = 0.0, tf_bott_cp = 0.0):
        
        self.bf_top = bf_top
        self.tf_top = tf_top
        self.t_web = t_web
        self.d_web = d_web
        self.bf_bott = bf_bott
        self.tf_bott = tf_bott
        
        self.bf_top_cp = bf_top_cp
        self.tf_top_cp = tf_top_cp
        self.bf_bott_cp = bf_bott_cp
        self.tf_bott_cp = tf_bott_cp
        
        self.pts = []
        self.calc_points()
        
        SP.__init__(self, self.pts)
        
        
    @property
    def bf_bott(self):
        return self.__bf_bott
    
    @bf_bott.setter
    def bf_bott(self, bf_bott):
        if bf_bott == None:
            self.__bf_bott = self.bf_top
        else:
            self.__bf_bott = bf_bott

    @property
    def tf_bott(self):
        return self.__tf_bott
    
    @tf_bott.setter
    def tf_bott(self, tf_bott):
        if tf_bott == None:
            self.__tf_bott = self.tf_top
        else:
            self.__tf_bott = tf_bott
        
        
    def calc_points(self):
        
        if self.bf_top_cp >= self.bf_bott_cp:
            
            d1 = self.bf_top_cp/2. - self.bf_bott_cp/2.
            d2 = self.bf_top_cp/2. - self.bf_bott/2.
            d3 = self.bf_top_cp/2. - self.bf_bott/2. + self.t_web
            d4 = self.tf_bott_cp + self.tf_bott
            d5 = self.bf_top_cp/2. + self.bf_bott/2. - self.t_web
            d6 = d4 + self.d_web
            d7 = d6 + self.tf_top
            d8 = d7 + self.tf_top_cp
            d9 = self.bf_top_cp/2. + self.bf_bott/2.
            d10 = self.bf_top_cp/2. + self.bf_bott_cp/2.
            
            self.pts.append((d1, 0.))
            self.pts.append((d1, self.tf_bott_cp))
            self.pts.append((d2, self.tf_bott_cp))
            self.pts.append((d3, d4))
            self.pts.append((d5, d4))
            self.pts.append((d5, d6))
            self.pts.append((d3, d6))
            self.pts.append((d3, d4))
            self.pts.append((d2, self.tf_bott_cp))
            self.pts.append((d2, d7))
            self.pts.append((0., d7))
            self.pts.append((0., d8))
            self.pts.append((self.bf_top_cp, d8))
            self.pts.append((self.bf_top_cp, d7))
            self.pts.append((d9, d7))
            self.pts.append((d9, self.tf_bott_cp))
            self.pts.append((d10, self.tf_bott_cp))
            
        elif self.bf_top_cp < self.bf_bott_cp:
            print("Not defined.")
            
    def __str__(self):
        """Provides summary of all properties."""
        
        user_input = (f'bf_top_cp = {self.bf_top_cp:.2f}, ' +
                      f'tf_top_cp = {self.tf_top_cp:.2f}, ' +
                      f'bf_top = {self.bf_top:.2f}, ' +
                      f'tf_top = {self.tf_top:.2f}, ' +
                      f't_web = {self.t_web:.2f}, ' +
                      f'd_web = {self.d_web:.2f}, ' +
                      f'bf_bott = {self.bf_bott:.2f}, ' +
                      f'tf_bott = {self.tf_bott:.2f}, ' +
                      f'bf_bott_cp = {self.bf_bott_cp:.2f}, ' +
                      f'tf_bott_cp = {self.tf_bott_cp:.2f} \n\n')
        
        all_props = super().__str__()
        
        return user_input + all_props
        
class ISP_AREMA(ISP):
    
    def __init__(self, bf_top, tf_top, t_web, d_web, 
                 bf_bott=None, tf_bott=None,
                 bf_top_cp=0.0, tf_top_cp=0.0, 
                 bf_bott_cp=0.0, tf_bott_cp=0.0):
        
        self.bf_top = bf_top
        self.tf_top = tf_top
        self.t_web = t_web
        self.d_web = d_web
        self.bf_bott = bf_bott
        self.tf_bott = tf_bott
        
        self.bf_top_cp = bf_top_cp
        self.tf_top_cp = tf_top_cp
        self.bf_bott_cp = bf_bott_cp
        self.tf_bott_cp = tf_bott_cp
        
        ISP.__init__(self, self.bf_top, self.tf_top, self.t_web, self.d_web,
                     self.bf_bott, self.tf_bott, 
                     self.bf_top_cp, self.tf_top_cp,
                     self.bf_bott_cp, self.tf_bott_cp)
        
        self.tsp_top = TSP(self.bf_top, self.tf_top, self.t_web, self.d_web,
                       self.bf_top_cp, self.tf_top_cp)
        self.tf_pl_top = RectSP(self.bf_top, self.tf_top, 
                         self.bf_top_cp, self.tf_top_cp)
        
        self.tsp_bott = TSP(self.bf_bott, self.tf_bott, self.t_web, self.d_web,
                       self.bf_bott_cp, self.tf_bott_cp)
        self.tf_pl_bott = RectSP(self.bf_bott, self.tf_bott, 
                         self.bf_bott_cp, self.tf_bott_cp)
        
        self.ryc_top = self.tsp_top.ry
        self.afc_top = self.tf_pl_top.area
        self.ryc_bott = self.tsp_bott.ry
        self.afc_bott = self.tf_pl_bott.area
        self.aw = self.t_web * self.d_web
        
#    covered by ISP parent class
#    def __str__(self):
#        """Provides summary of all properties."""
        
    
class Box_AREMA(Box):
    
    def __init__(self, bf_top, tf_top, t_web, d_web, bf_bott = None, 
                 tf_bott = None, bf_top_cp = 0.0, tf_top_cp = 0.0,
                 bf_bott_cp = 0.0, tf_bott_cp = 0.0):
        
        self.bf_top = bf_top
        self.tf_top = tf_top
        self.t_web = t_web
        self.d_web = d_web
        self.bf_bott = bf_bott
        self.tf_bott = tf_bott
        
        self.bf_top_cp = bf_top_cp
        self.tf_top_cp = tf_top_cp
        self.bf_bott_cp = bf_bott_cp
        self.tf_bott_cp = tf_bott_cp
        
        Box.__init__(self, self.bf_top, self.tf_top, self.t_web, self.d_web,
             self.bf_bott, self.tf_bott, self.bf_top_cp, self.tf_top_cp,
             self.bf_bott_cp, self.tf_bott_cp)
        
        self.sum_st = self.calc_sum_st()
        self.enclosed_area = self.calc_enclosed_area()
        self.aw = 2*self.t_web*self.d_web
        
    def calc_sum_st(self):
        bf_top_avg = (self.bf_top_cp + self.bf_top)/2
        tf_top = self.tf_top_cp + self.tf_top
        bf_bott_avg = (self.bf_bott_cp + self.bf_bott)/2
        tf_bott = self.tf_bott_cp + self.tf_bott
        
        r1 = bf_top_avg/tf_top
        r2 = self.d_web/self.t_web
        r3 = bf_bott_avg/tf_bott
        
        return r1 + 2*r2 + r3
    
    def calc_enclosed_area(self):
        #assumes the box is rectangular
        top = self.bf_top - 2*self.t_web
        side = self.d_web + self.tf_top/2. + self.tf_bott/2.
        
        return top*side