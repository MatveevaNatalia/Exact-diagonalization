import unittest
from Func import *

class TestGround_State(unittest.TestCase):
    def test_6_5(self):
        result = Ground_State(6,5)
        self.assertEqual(result, [1, 1, 1, 1, 1, 0])

    def test_6_6(self):
        result = Ground_State(6,6)
        self.assertEqual(result, [1, 1, 1, 1, 1, 1])
        
    def test_error_smaller(self):
        with self.assertRaises(ValueError):
            Ground_State(5,6)
            
    def test_error_negative(self):
        with self.assertRaises(ValueError):
            Ground_State(-5,2)        
        with self.assertRaises(ValueError):
            Ground_State(5,-2)        

class TestTo_integer(unittest.TestCase):
    def test_7(self):
        result = To_integer([1,1,1])
        self.assertEqual(result, 7)

class TestKillState(unittest.TestCase):
    def test_4(self):
        res_state, res_coeff = Kill_State(0, 12, 4)
        self.assertEqual(res_state, 4)
        
    def test_too_big_index(self):
         with self.assertRaises(ValueError):
            Kill_State(5,12,4) 
            
    def test_negative_numbers(self):
        with self.assertRaises(ValueError):
            Kill_State(-1,12,4)              
        with self.assertRaises(ValueError):
            Kill_State(1,-12,4) 
        with self.assertRaises(ValueError):
            Kill_State(1,12,-4)    
        with self.assertRaises(ValueError):
            Kill_State(-1,12,-4)     

class TestCreateState(unittest.TestCase):
    def test_12(self):
        res_state, res_coeff = Create_State(0, 4, 4)
        self.assertEqual(res_state, 12)
        
    def test_too_big_index(self):
        with self.assertRaises(ValueError):
            Create_State(5,12,4) 
            
    def test_negative_numbers(self):
        with self.assertRaises(ValueError):
            Create_State(-1,12,4)              
        with self.assertRaises(ValueError):
            Create_State(1,-12,4) 
        with self.assertRaises(ValueError):
            Create_State(1,12,-4)    
        with self.assertRaises(ValueError):
            Create_State(-1,12,-4) 
            
if __name__ == '__main__':
    unittest.main()