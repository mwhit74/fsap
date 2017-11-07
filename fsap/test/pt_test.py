import unittest
from utils import point as pt

class ThreeDPointTest(unittest.TestCase):
    def setUp(self):
        self.three_d_pt = pt.ThreeDPoint(2.3, 4.5, 9.4)

    def tearDown(self):
        del(self.three_d_pt)

    def test_return_x(self):
        self.assertEqual(self.three_d_pt.x, 2.3, 'Incorrect x coord')

    def test_return_y(self):
        self.assertEqual(self.three_d_pt.y, 4.5, 'Incorrect y coord')

    def test_return_z(self):
        self.assertEqual(self.three_d_pt.z, 9.4, 'Incorrect z coord')

    def test_string(self):
        self.assertEqual(str(self.three_d_pt), '2.300 4.500 9.400', 'Incorrect string point')

    def test_return_rect_3D(self):
        self.assertEqual(self.three_d_pt.rect_3D(), (2.3, 4.5, 9.4), 'Incorrect tuple')

    def test_return_rect_2D(self):
        self.assertEqual(self.three_d_pt.rect_2D(), (2.3, 4.5), 'Incorrect tuple')

    def test_set_new_x(self):
        self.three_d_pt.x = 3.7
        self.assertEqual(self.three_d_pt.x, 3.7, 'Incorrect new x')

    def test_set_new_y(self):
        self.three_d_pt.y = 5.9
        self.assertEqual(self.three_d_pt.y, 5.9, 'Incorrect new y')

    def test_set_new_z(self):
        self.three_d_pt.z = 8.1
        self.assertEqual(self.three_d_pt.z, 8.1, 'Incorrect new z')

