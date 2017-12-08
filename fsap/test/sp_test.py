import unittest
from fsap.section_property import section_property as sp

class SectionPropertyTest(unittest.TestCase):
    def setUp(self):
        self.points1 = [(0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)]
        self.points2 = [(0.0,0.0),(1.0,1.0),(1.0,0.0),(0.0,1.0)]


    def test_centroid(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.centroid.x, 0.5, "Incorrect center x-coord")
        self.assertEqual(sp1.centroid.y, 0.5, "Incorrect center y-coord")


    def test_order_points(self):
        sp2 = sp.SectionProperty(self.points2)
        self.assertItemsEqual(sp2.points[0].cc(),(0.0, 1.0), "Incorrect first point")
        self.assertItemsEqual(sp2.points[1].cc(),(0.0, 0.0), "Incorrect second point")
        self.assertItemsEqual(sp2.points[2].cc(),(1.0, 0.0), "Incorrect third point")
        self.assertItemsEqual(sp2.points[3].cc(),(1.0, 1.0), "Incorrect fourth point")
        self.assertItemsEqual(sp2.points[4].cc(),(0.0, 1.0), "Incorrect fifth point")


    def test_max_y(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.max_y, 1.0, "Incorrect max y")


    def test_max_x(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.max_x, 1.0, "Incorrect max x " +
                str(sp1.max_x))


    def test_min_y(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.min_y, 0.0, "Incorrect min y")


    def test_min_x(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.min_x, 0.0, "Incorrect min x")


    def test_bounding_box(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertItemsEqual(sp1.box,(1.0, 1.0, 0.0, 0.0), "Incorrect box")


    def test_width(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.width, 1.0, "Incorrect width")


    def test_height(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.height, 1.0, "Incorrect height")


    def test_area(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.area, 1.0, "Incorrect area " + str(sp1.area))


    def test_ena_x(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.ena_x, 0.5, "Incorrect ena_x")


    def test_ena_y(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.ena_y, 0.5, "Incorrect ena_y")


    def test_ixx_x(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.ixx_x, 0.333, "Incorrect ixx_x " + str(sp1.ixx_x))


    def test_iyy_y(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.iyy_y, 0.333, "Incorrect iyy_y " + str(sp1.iyy_y))

        
    def test_ixx_c(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.ixx_c, 0.083, "Incorrect ixx_c " + str(sp1.ixx_c))


    def test_iyy_c(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.iyy_c, 0.083, "Incorrect iyy_c " + str(sp1.iyy_c))


    def test_ixy_xy(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.ixy_xy, 0.25, "Incorrect ixy_xy " +
                str(sp1.ixy_xy))


    def test_ixy_c(self):
        sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(sp1.ixy_c, 0.0, "Incorrect ixy_x " + str(sp1.ixy_c))
