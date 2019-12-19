import unittest
from fsap.anen import model as ml

class ModelTest(unittest.TestCase):
    def setUp(self):
        sup = [(1,1,1),
               (3,0,1),
               (4,0,1)]

    def test_num_dof(self):
        
