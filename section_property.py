"""Section Properties"""

def area(xy_coords):
	return -1*loop(xy_coords, area_eq)
	
def area_eq(cur_x, cur_y, next_x, next_y):
	return (next_y - cur_y)*(next_x + cur_x)/2.0

	
def ENA_x(xy_coords):
	A = area(xy_coords)
	return -1.0/A*loop(xy_coords, ENA_x_eq)
	
def ENA_x_eq(cur_x, cur_y, next_x, next_y):
	return ((next_y - cur_y)/8.0)*((next_x + cur_x)**2.0 + (next_x - cur_x)**2.0/3.0)
	

def ENA_y(xy_coords):
	A = area(xy_coords)
	return 1.0/A*loop(xy_coords, ENA_y_eq)
	
def ENA_y_eq(cur_x, cur_y, next_x, next_y):
	return ((next_x - cur_x)/8.0)*((next_y + cur_y)**2.0 + (next_y - cur_y)**2.0/3.0)
	
	
def I_about_x(xy_coords):
	return loop(xy_coords, I_about_x_eq)
	
def I_about_x_eq(cur_x, cur_y, next_x, next_y):
	return ((next_x - cur_x)*(next_y + cur_y)/24.0)*((next_y + cur_y)**2.0 + (next_y - cur_y)**2.0)
	
	
def I_about_y(xy_coords):
	return -1*loop(xy_coords, I_about_y_eq)
	
def I_about_y_eq(cur_x, cur_y, next_x, next_y):
	return ((next_y - cur_y)*(next_x + cur_x)/24.0)*((next_x + cur_x)**2.0 + (next_x - cur_x)**2.0)
	
	
def I_about_x_centroid(xy_coords):
	A = area(xy_coords)
	y_bar = ENA_y(xy_coords)
	I_x = I_about_x(xy_coords)
	
	return (I_x - A*y_bar**2)
	

def I_about_y_centroid(xy_coords):
	A = area(xy_coords)
	x_bar = ENA_x(xy_coords)
	I_y = I_about_y(xy_coords)
	
	return (I_y - A*x_bar**2)

	
def loop(xy_coords, func):
	var = 0
	
	cur_x = xy_coords[0][0]
	cur_y = xy_coords[0][1]
	
	for xy in xy_coords:
		next_x = xy[0]
		next_y = xy[1]
		
		var = var + func(cur_x, cur_y, next_x, next_y)
			
		cur_x = next_x
		cur_y = next_y
	
	return var
	
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