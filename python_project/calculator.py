"""
A simple calculator module for demonstration purposes.
"""

def add(a, b):
    """
    Add two numbers together.
    
    Args:
        a (int or float): The first number.
        b (int or float): The second number.
        
    Returns:
        int or float: The sum of a and b.
    """
    return a + b

def divide(a, b):
    """
    Divide the first number by the second number.
    
    Args:
        a (int or float): The dividend.
        b (int or float): The divisor. Must not be zero.
        
    Returns:
        float: The result of the division.
        
    Raises:
        ValueError: If b is zero.
    """
    if b == 0:
        raise ValueError("Cannot divide by zero")
    return a / b