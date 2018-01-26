import unittest
from fsap.geom import point as pt

class Point2DTest(unittest.TestCase):
    def setUp(self):
        self.two_d_pt = pt.Point2D(2.3, 4.5)

    def tearDown(self):
        del(self.two_d_pt)

    def test_return_x(self):
        self.assertEqual(self.two_d_pt.x, 2.3, 'Incorrect x coord')

    def test_return_y(self):
        self.assertEqual(self.two_d_pt.y, 4.5, 'Incorrect y coord')

    def test_string(self):
        self.assertEqual(str(self.two_d_pt), '(2.300,4.500)', 'Incorrect string point')

    def test_coord(self):
        self.assertEqual(self.two_d_pt.cc()[0], 2.3, 'Incorrect tuple')
        self.assertEqual(self.two_d_pt.cc()[1], 4.5, 'Incorrect tuple')

    def test_set_new_x(self):
        self.two_d_pt.x = 3.7
        self.assertEqual(self.two_d_pt.x, 3.7, 'Incorrect new x')

    def test_set_new_y(self):
        self.two_d_pt.y = 5.9
        self.assertEqual(self.two_d_pt.y, 5.9, 'Incorrect new y')
    

class Point3DTest(unittest.TestCase):
    def setUp(self):
        self.three_d_pt = pt.Point3D(2.3, 4.5, 9.4)

    def tearDown(self):
        del(self.three_d_pt)

    def test_return_x(self):
        self.assertEqual(self.three_d_pt.x, 2.3, 'Incorrect x coord')

    def test_return_y(self):
        self.assertEqual(self.three_d_pt.y, 4.5, 'Incorrect y coord')

    def test_return_z(self):
        self.assertEqual(self.three_d_pt.z, 9.4, 'Incorrect z coord')

    def test_string(self):
        self.assertEqual(str(self.three_d_pt), '(2.300,4.500,9.400)', 'Incorrect string point')

    def test_coord(self):
        self.assertEqual(self.three_d_pt.cc()[0], 2.3, 'Incorrect tuple')
        self.assertEqual(self.three_d_pt.cc()[1], 4.5, 'Incorrect tuple')
        self.assertEqual(self.three_d_pt.cc()[2], 9.4, 'Incorrect tuple')

    def test_set_new_x(self):
        self.three_d_pt.x = 3.7
        self.assertEqual(self.three_d_pt.x, 3.7, 'Incorrect new x')

    def test_set_new_y(self):
        self.three_d_pt.y = 5.9
        self.assertEqual(self.three_d_pt.y, 5.9, 'Incorrect new y')

    def test_set_new_z(self):
        self.three_d_pt.z = 8.1
        self.assertEqual(self.three_d_pt.z, 8.1, 'Incorrect new z')


