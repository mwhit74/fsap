import math
from fsap.utils.point import Point2D
import pdb


def less(a, b, center):
    """Comparison function to determine the order of two points"""
    #if a is left of center and b is right of center, a is ccw from b
    if a.x - center.x >= 0.0 and b.x - center.x < 0.0:
        return True
    #if a is right of center and be is left of center, a is cw from b
    if a.x - center.x < 0.0 and b.x - center.x >= 0:
        return False
    #if a, b, and center lie on the same same vertical line
    if a.x - center.x == 0.0 and b.x - center.x == 0.0:
        #if a.y is greater than b.y, a.y is ccw from b.y otherwise a.y is cw
        #from b.y
        return a.y > b.y

        #if not requiring points to be in first quadrant
        #if a.y - center.y >= 0 or b.y - center.y >=0:
        #    return a.y > b.y
        #return b.y > a.y

    #if a.x and b.x are on the same side of center calculate the
    #cross-product  of (center -> a) x (center -> b)
    det = ((a.x - center.x)*(b.y - center.y) - 
            (b.x - center.x)*(a.y - center.y))
    #if the cross-product is positive a is ccw from b, otherwise a is cw
    #from b
    if (det < 0):
        return True
    else:
        return False

    #a and b lie on the same line from the center
    #if a is farther than b, a is ccw from b, otherwise a is cw from b
    d1 = math.pow((a.x - center.x),2) + math.pow((a.y - center.y),2)
    d2 = math.pow((b.x - center.x),2) + math.pow((b.y - center.y),2)
    return d1 > d2

if __name__ == "__main__":
    points1 = [Point2D(0.0,0.0),Point2D(1.0,0.0),Point2D(1.0,1.0),Point2D(0.0,1.0),Point2D(0.0, 0.0)]
    points2 = [Point2D(0.0,0.0),Point2D(1.0,1.0),Point2D(1.0,0.0),Point2D(0.0,1.0),Point2D(0.0, 0.0)]

    #pdb.set_trace()
    print less(points1[0], points1[1], Point2D(0.5, 0.5))
    print less(points2[1], points2[2], Point2D(0.5, 0.5))
