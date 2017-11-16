import unittest
from fsap.section_property import section_property as sp

class SectionPropertyTest(unittest.TestCase):
    def setUp(self):
        self.sp1 = sp.SectionProperty()
        self.points = [(0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0), (0.0, 0.0)]

    def test_convert_points(self):
        self.sp1.convert_to_points(self.points)
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
