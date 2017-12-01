import unittest
from fsap.section_property import section_property as sp

class SectionPropertyTest(unittest.TestCase):
    def setUp(self):
        self.points1 = [(0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0), (0.0, 0.0)]
        self.points2 = [(0.0,0.0),(1.0,1.0),(1.0,0.0),(0.0,1.0), (0.0, 0.0)]


    def test_convert_points(self):
        self.sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(self.sp1.points[0].x, 0.0, "Incorrect x-coord point 1")
        self.assertEqual(self.sp1.points[0].y, 0.0, "Incorrect y-coord point 1")
        self.assertEqual(self.sp1.points[1].x, 1.0, "Incorrect x-coord point 2")
        self.assertEqual(self.sp1.points[1].y, 0.0, "Incorrect y-coord point 2")
        self.assertEqual(self.sp1.points[2].x, 1.0, "Incorrect x-coord point 3")
        self.assertEqual(self.sp1.points[2].y, 1.0, "Incorrect y-coord point 3")
        self.assertEqual(self.sp1.points[3].x, 0.0, "Incorrect x-coord point 4")
        self.assertEqual(self.sp1.points[3].y, 1.0, "Incorrect y-coord point 4")
        self.assertEqual(self.sp1.points[4].x, 0.0, "Incorrect x-coord point 5")
        self.assertEqual(self.sp1.points[4].y, 0.0, "Incorrect y-coord point 5")

    def test_center(self):
        self.sp1 = sp.SectionProperty(self.points1)
        self.assertEqual(self.sp1.center.x, 0.5, "Incorrect center x-coord")
        self.assertEqual(self.sp1.center.y, 0.5, "Incorrect center y-coord")

    def test_order_points(self):
        self.sp1 = sp.SectionProperty(self.points2)
        self.assertItemsEqual(self.sp1.points[0].cc(),(0.0, 0.0), "Incorrect first point")
        self.assertItemsEqual(self.sp1.points[1].cc(),(1.0, 0.0), "Incorrect second point")
        self.assertItemsEqual(self.sp1.points[2].cc(),(1.0, 1.0), "Incorrect third point")
        self.assertItemsEqual(self.sp1.points[3].cc(),(0.0, 1.0), "Incorrect fourth point")
        self.assertItemsEqual(self.sp1.points[4].cc(),(0.0, 0.0), "Incorrect fifth point")

