"""
Unit tests for the calculator module.
"""

import unittest
from calculator import add, divide

class TestCalculator(unittest.TestCase):
    
    def test_add_positive_numbers(self):
        """Test adding two positive numbers."""
        self.assertEqual(add(2, 3), 5)
        
    def test_add_negative_numbers(self):
        """Test adding two negative numbers."""
        self.assertEqual(add(-2, -3), -5)
        
    def test_add_positive_and_negative(self):
        """Test adding a positive and a negative number."""
        self.assertEqual(add(5, -3), 2)
        
    def test_add_zero(self):
        """Test adding zero to a number."""
        self.assertEqual(add(0, 5), 5)
        self.assertEqual(add(5, 0), 5)
        
    def test_divide_positive_numbers(self):
        """Test dividing two positive numbers."""
        self.assertAlmostEqual(divide(10, 3), 3.333333333333333)
        
    def test_divide_by_zero_raises_exception(self):
        """Test that dividing by zero raises a ValueError."""
        with self.assertRaises(ValueError):
            divide(10, 0)
            
    def test_divide_negative_numbers(self):
        """Test dividing with negative numbers."""
        self.assertEqual(divide(-10, 2), -5)
        self.assertEqual(divide(10, -2), -5)
        self.assertEqual(divide(-10, -2), 5)

if __name__ == '__main__':
    unittest.main()