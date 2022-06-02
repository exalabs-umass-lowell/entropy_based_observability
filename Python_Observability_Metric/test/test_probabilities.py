import unittest
import sys

from random import Random
random = Random(1337)

sys.path.append('.')
sys.path.append('./src')

# local files
from src.probabilities import *

class TestSetup(unittest.TestCase):

    def test_calculate_prob_ff_from_recent_motion_should_return_min(self):
        len_of_time_window = 16
        recent_motion_data = np.zeros(len_of_time_window).astype(float)
        res = calculate_prob_ff_from_recent_motion(recent_motion_data)
        self.assertEqual(res, 0.0)

    def test_calculate_prob_ff_from_recent_motion_should_return_max(self):
        len_of_time_window = 16
        recent_motion_data = np.ones(len_of_time_window).astype(float)
        res = calculate_prob_ff_from_recent_motion(recent_motion_data)
        self.assertEqual(res, 1.0)
        
    def test_calculate_prob_ff_from_recent_motion_should_return_in_between(self):
        len_of_time_window = 16
        half_len = int(np.floor(len_of_time_window/2))
        one_half = np.ones(half_len).astype(float)
        zero_half = np.zeros(half_len).astype(float)
        recent_motion_data = np.random.permutation(np.append(one_half,zero_half))
        res = calculate_prob_ff_from_recent_motion(recent_motion_data)
        self.assertGreater(res, 0.0)
        self.assertLess(res, 1.0)

    def test_calculate_prob_dist_agents_in_sites(self):
        # for full
        influential_range = 4
        num_local_agents = 4
        res = calculate_prob_dist_agents_in_sites(num_local_agents, influential_range)
        self.assertEqual(len(res), 1)
        # for empty
        influential_range = 4
        num_local_agents = 0
        res = calculate_prob_dist_agents_in_sites(num_local_agents, influential_range)
        self.assertEqual(len(res), 1)
        
    def test_calculate_prob_dist_agents_out_of_range(self):
        num_agents_out_of_range = 2
        num_sites_out_of_range = 60
        res = calculate_prob_dist_agents_out_of_range(num_agents_out_of_range, num_sites_out_of_range)
        self.assertEqual(res.shape[0], num_agents_out_of_range)
        
        
    # def test_calculate_prob_dist_agents_out_of_range_should_equal_agents_in_sites(self):
    #     num_agents_out_of_range = 4
    #     num_sites_out_of_range = 20
    #     res = calculate_prob_dist_agents_out_of_range(num_agents_out_of_range, num_sites_out_of_range)
    #     res2 = calculate_prob_dist_agents_in_sites(num_agents_out_of_range, num_sites_out_of_range)
    #     self.assertEqual(res, res2)
        
    def test_calculate_agents_energy_level_count(self):
        num_agents_out_of_range = 1
        num_sites_out_of_range = 60
        res = calculate_agents_energy_level_count(num_agents_out_of_range, num_sites_out_of_range)
        self.assertEqual(res.shape[0], num_agents_out_of_range)
        

    def test_calculate_energy_distribution_other_sites(self):
        num_of_agents_out_of_range = random.randint(0, 20)
        num_of_sites_out_of_range = random.randint(20, 40)
        prob_dist = calculate_energy_distribution_other_sites(num_of_agents_out_of_range, num_of_sites_out_of_range)
        self.assertEqual(prob_dist.shape[0], num_of_agents_out_of_range)
        
if __name__ == '__main__':
    unittest.main()