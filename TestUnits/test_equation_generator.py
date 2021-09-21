""" Contains the tests for the EquationGenerator class.
"""

# Set the current working directory to the source diretory.
import os
import copy as cp
path1 = cp.deepcopy(os.getcwd())
os.chdir(os.path.dirname(__file__) + os.sep + "..")

# Imports: The unittest module.
import unittest

# Imports: Class to be tested.
from Program.mathematica_generator import EquationGenerator


class EquationGenerator(unittest.TestCase):

    def test_sum(self):
        self.assertEqual(sum([1, 2, 3]), 6, "Should be 6")


if __name__ == '__main__':
    # Run the tests.
    unittest.main()

    # Set the working directory back to the original.
    os.chdir(path1)
