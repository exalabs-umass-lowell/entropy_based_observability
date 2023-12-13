import unittest
import sys

from random import Random
random = Random(1337)

sys.path.append('.')
sys.path.append('./src')

# local files
from src.simulation import *

class TestSetup(unittest.TestCase):

    def test_step_configuration(self):
        return True

    def test_modify_start_pose(self):
        return True

    def test_state_prediction_from_conditional_prob(self):
        return True

    def test_perform_simulation(self):
        return True

    def test_statistical_mechanics_based_traffic(self):      
        return True

if __name__ == '__main__':
    unittest.main()