from Points import Point2D

class SectionProperty:
    def __init__(xy_coords):
        self.points = convert_to_points(xy_coords)

        self.area = self.area()
        self.ena_x = self.ENA_x()
        self.ena_y = self.ENA_y()
        self.ixx_x = self.I_about_x()
        self.iyy_y = self.I_about_y()
        self.ixx_c = self.I_about_x_centroid()
        self.iyy_c = self.I_about_y_centroid()

    def convert_to_points(xy_coords):
        points = []
        for xy in xy_coords:
            points.append(Point2D(xy[0], xy[1]))
        return points

    def area():
            return -1*loop(self.points, self.area_eq)
            
    def area_eq(cur_x, cur_y, next_x, next_y):
            return (next_y - cur_y)*(next_x + cur_x)/2.0

            
    def ENA_x():
            A = self.area
            return -1.0/A*loop(self.points, self.ENA_x_eq)
            
    def ENA_x_eq(cur_x, cur_y, next_x, next_y):
            return ((next_y - cur_y)/8.0)*((next_x + cur_x)**2.0 + (next_x - cur_x)**2.0/3.0)
            

    def ENA_y():
            A = self.area
            return 1.0/A*loop(self.points, self.ENA_y_eq)
            
    def ENA_y_eq(cur_x, cur_y, next_x, next_y):
            return ((next_x - cur_x)/8.0)*((next_y + cur_y)**2.0 + (next_y - cur_y)**2.0/3.0)
            
            
    def I_about_x():
            return loop(self.points, self.I_about_x_eq)
            
    def I_about_x_eq(cur_x, cur_y, next_x, next_y):
            return ((next_x - cur_x)*(next_y + cur_y)/24.0)*((next_y + cur_y)**2.0 + (next_y - cur_y)**2.0)
            
            
    def I_about_y():
            return -1*loop(self.points, self.I_about_y_eq)
            
    def I_about_y_eq(cur_x, cur_y, next_x, next_y):
            return ((next_y - cur_y)*(next_x + cur_x)/24.0)*((next_x + cur_x)**2.0 + (next_x - cur_x)**2.0)
            
            
    def I_about_x_centroid():
            A = self.area
            y_bar = self.ENA_y
            I_x = self.I_about_x()
            
            return (I_x - A*y_bar**2)
            

    def I_about_y_centroid(xy_coords):
            A = self.area
            x_bar = self.ENA_x
            I_y = self.I_about_y
            
            return (I_y - A*x_bar**2)

            
    def loop(self.points, func):
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
	print "ENA in x: " + str(ENA_x(xy_coords))
	print "ENA in y: " + str(ENA_y(xy_coords))
	print "I about x: " + str(I_about_x(xy_coords))
	print "I about y: " + str(I_about_y(xy_coords))
	print "I about x at centroid: " + str(I_about_x_centroid(xy_coords))
	print "I about y at centroid: " + str(I_about_y_centroid(xy_coords))
