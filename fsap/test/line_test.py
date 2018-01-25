import unittest
from fsap.geom import point as pt
from fsap.geom import line as ln

class Line2DTest(unittest.TestCase):
    def setUp(self):
        self.pt1 = pt.Point2D(2.3, 4.5)
        self.pt2 = pt.Point2D(7.8, 12.6)
        self.ln2d = ln.Line2D(self.pt1, self.pt2)

    def test_return_x1(self):
        self.assertEqual(self.ln2d.pt1.x, 2.3, "Incorrect x1 coord")

    def test_return_y1(self):
        self.assertEqual(self.ln2d.pt1.y, 4.5, "Incorrect y1 coord")

    def test_return_x2(self):
        self.assertEqual(self.ln2d.pt2.x, 7.8, "Incorrect x2 coord")

    def test_return_y2(self):
        self.assertEqual(self.ln2d.pt2.y, 12.6, "Incoorect y2 coord")

    def test_str(self):
        self.assertEqual(str(self.ln2d), "(2.3,4.5) (7.8,12.6)", "Incorrect str")

    def test_length(self):
        self.assertAlmostEqual(self.ln2d.length(), 9.7908, places=4, msg="Incorrect length")

class Line3DTest(unittest.TestCase):
    def setUp(self):
        self.pt1 = pt.Point3D(2.3, 4.5, 1.5)
        self.pt2 = pt.Point3D(7.8, 12.6, 15.9)
        self.ln3d = ln.Line3D(self.pt1, self.pt2)

    def test_return_x1(self):
        self.assertEqual(self.ln3d.pt1.x, 2.3, "Incorrect x1 coord")

    def test_return_y1(self):
        self.assertEqual(self.ln3d.pt1.y, 4.5, "Incorrect y1 coord")

    def test_return_z1(self):
        self.assertEqual(self.ln3d.pt1.z, 1.5, "Incorrect z1 coord")

    def test_return_x2(self):
        self.assertEqual(self.ln3d.pt2.x, 7.8, "Incorrect x2 coord")

    def test_return_y2(self):
        self.assertEqual(self.ln3d.pt2.y, 12.6, "Incoorect y2 coord")

    def test_return_z2(self):
        self.assertEqual(self.ln3d.pt2.z, 15.9, "Incorrect z2 coord")

    def test_str(self):
        self.assertEqual(str(self.ln3d), "(2.3,4.5,1.5) (7.8,12.6,15.9)", "Incorrect str")

    def test_length(self):
        self.assertAlmostEqual(self.ln3d.length(), 17.4132, places=4, msg="Incorrect length")
