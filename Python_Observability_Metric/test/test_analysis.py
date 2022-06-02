import unittest
import sys

from random import Random
random = Random(1337)

sys.path.append('.')
sys.path.append('./src')

# local files
from src.utility import *
from src.analysis import *

# Traffic Model prameters
# and Parameters for information quantification
global_params = {
    "num_of_sites": 64,
    "num_of_agents": 20,
    "interaction_coeff": 4.0,
    "external_field_coeff": 1.5,
    "len_of_time": 64,
    "influential_range": 4,
    "len_of_time_window": 16
}
class TestSetup(unittest.TestCase):

    def test_state_prediction_from_neighboring_site_states(self):
        return True

    def test_derive_and_normalize_likelihood_from_bayesian(self):
        return True

    def test_motion_data_filtering(self):
        return True

    # def test_calculate_probs_out_of_range_agents(self):
    #     global_params = {
    #         "num_of_sites": 64,
    #         "num_of_agents": 5,
    #         "influential_range": 4,
    #     }
    #     num_of_local_agents = global_params["num_of_agents"] - 1
    #     num_of_out_of_range_agents = global_params["num_of_agents"] - num_of_local_agents
    #     num_of_out_of_range_sites = global_params["num_of_sites"] - global_params["influential_range"]
    #     res = calculate_prob_dist_agents_out_of_range(num_of_out_of_range_agents, num_of_out_of_range_sites)
    #     self.assertEqual(len(res), 0)

    # def test_calculate_energy_distribution_other_sites(self):
    #     # given
    #     num_agents = 19;
    #     num_sites = 59; 
    #     # when
    #     energy_dist = calculate_energy_distribution_other_sites(num_agents, num_sites)
    #     # then
    #     self.assertEqual(len(energy_dist), num_agents)

    
    def test_calculate_direct_global_entropy_is_same_as_estimate(self):
        # given simulation parameters for a small enough system
        global_params = {
            "num_of_sites": 20,
            "num_of_agents": 5,
            "influential_range": 4,
        }
        # when the entropy of the global system is calculated using either method
        entropy1 = estimate_global_entropy(global_params)
        entropy2 = calculate_direct_global_entropy(global_params)
        # then the two methods should return the same value
        self.assertEqual(entropy1, entropy2)
        
        
    def test_calculate_direct_global_entropy_is_zero_for_full_case(self):
        # given simulation parameters for a small enough system
        global_params = {
            "num_of_sites": 20,
            "num_of_agents": 20,
            "influential_range": 4,
        }
        # when the entropy is calculated directly for a full system
        entropy = calculate_direct_global_entropy(global_params)
        # then it should be 0.0
        self.assertEqual(entropy, 0.0)
        
    def test_calculate_direct_global_entropy_is_zero_for_empty_case(self):
        # given simulation parameters for a small enough system
        global_params = {
            "num_of_sites": 20,
            "num_of_agents": 0,
            "influential_range": 4,
        }
        # when the entropy is calculated directly
        entropy = calculate_direct_global_entropy(global_params)
        # then it should be 0.0
        self.assertEqual(entropy, 0.0)

    def test_calculate_direct_global_entropy_is_non_zero_for_intermediate_case(self):
        # given simulation parameters for a small enough system
        global_params = {
            "num_of_sites": 64,
            "num_of_agents": 10,
            "influential_range": 4,
        }
        # when the entropy is calculated directly
        entropy = calculate_direct_global_entropy(global_params)
        # it should be non-zero
        self.assertNotEqual(entropy, 0.0)
        
        
    # def test_estimate_global_entropy_is_zero_for_full_case(self):
    #     # given simulation parameters for a small enough system
    #     global_params = {
    #         "num_of_sites": 20,
    #         "num_of_agents": 20,
    #         "influential_range": 4,
    #     }
    #     # when the entropy is calculated directly for a full system
    #     entropy = estimate_global_entropy(global_params)
    #     # then it should be 0.0
    #     self.assertEqual(entropy, 0.0)
        
    # def test_estimate_global_entropy_is_zero_for_empty_case(self):
    #     # given simulation parameters for a small enough system
    #     global_params = {
    #         "num_of_sites": 20,
    #         "num_of_agents": 0,
    #         "influential_range": 4,
    #     }
    #     # when the entropy is calculated directly
    #     entropy = estimate_global_entropy(global_params)
    #     # then it should be 0.0
    #     self.assertEqual(entropy, 0.0)

    # def test_estimate_global_entropy_is_zero_for_empty_case(self):
    #     # given simulation parameters for a small enough system
    #     global_params = {
    #         "num_of_sites": 20,
    #         "num_of_agents": 0,
    #         "influential_range": 4,
    #     }
    #     # when the entropy is calculated directly
    #     entropy = estimate_global_entropy(global_params)
    #     # then it should be 0.0
    #     self.assertEqual(entropy, 0.0)

    # def test_estimate_global_entropy_is_non_zero_far_minimun_case(self):
    #     # given simulation parameters for a small enough system
    #     global_params = {
    #         "num_of_sites": 20,
    #         "num_of_agents": 10,
    #         "influential_range": 4,
    #     }
    #     # when the entropy is calculated directly
    #     entropy = estimate_global_entropy(global_params)
    #     # it should be non-zero
    #     self.assertNotEqual(entropy, 0.0)
        
    def test_calculate_probs_dist_global_sites_given_local_agents(self):
        num_local_agents = 0
        global_params = {
            "num_of_sites": 20,
            "num_of_agents": 10,
            "influential_range": 4
        }
        probs_dist = calculate_probs_dist_global_sites_given_local_agents(num_local_agents, global_params)
        print(probs_dist)
        self.assertEquals(len(probs_dist), 1)     
    
    def test_calculate_agents_energy_level_count(self):
        num_of_agents = 2
        num_of_sites = 6
        res = calculate_agents_energy_level_count(num_of_agents, num_of_sites)
        self.assertEqual(len(res), num_of_agents)

    def test_calculate_probs_ff_given_possible_local_sites(self):
        global_params = {
            "num_of_sites": 20,
            "num_of_agents": 5,
            "interaction_coeff": 4.0,
            "external_field_coeff": 1.5,
            "influential_range": 4
        }
        possible_local_sites = calculate_possible_sites(global_params["influential_range"])
        probs_ff = calculate_probs_ff_given_possible_local_sites(possible_local_sites, global_params)
        # print(probs_ff)
        # print(possible_local_sites)
        self.assertEqual(len(probs_ff), len(possible_local_sites))
        
    # def test_calculate_likelihood_local_sites_given_agent_motion(self):
    #     global_params = {
    #         "num_of_sites": 6,
    #         "num_of_agents": 3,
    #         "interaction_coeff": 4.0,
    #         "external_field_coeff": 1.5,
    #         "len_of_time": 64,
    #         "influential_range": 3,
    #         "len_of_time_window": 4
    #     }
    #     recent_agent_motion = np.zeros(global_params["len_of_time_window"]).astype(float)
    #     likeliehood_vec = calculate_likelihood_local_sites_given_agent_motion(recent_agent_motion, global_params)
    #     # then it produce a probability vec which should sum to 1
    #     self.assertEqual(np.sum(likeliehood_vec), 1.0)


    def test_calculate_prob_dist_remote_sites_given_local_sites(self):
        global_params = {
            "num_of_sites": 6,
            "num_of_agents": 3,
            "interaction_coeff": 4.0,
            "external_field_coeff": 1.5,
            "len_of_time": 64,
            "influential_range": 3,
            "len_of_time_window": 4
        }
        local_sites = np.zeros(3).astype(float)
        #likeliehood_vec = calculate_likelihood_local_sites_given_agent_motion(recent_agent_motion, global_params)
        possible_remote_sites, prob_dist = calculate_prob_dist_remote_sites_given_local_sites_and_agent_num(local_sites, global_params["num_of_agents"], global_params)
        # then it produce a probability vec which should sum to 1
        self.assertEqual(np.sum(prob_dist), 1.0)

        

    # def test_calculate_observability(self):
    #     global_params = {
    #         "num_of_sites": 64,
    #         "num_of_agents": 5,
    #         "interaction_coeff": 4.0,
    #         "external_field_coeff": 1.5,
    #         "len_of_time": 64,
    #         "influential_range": 4,
    #         "len_of_time_window": 16
    #     }
    #     len_of_time_window = global_params["len_of_time_window"]
    #     recent_motion_data = np.ones(len_of_time_window).astype(float)
    #     global_entropy = 12
    #     res = calculate_observability(recent_motion_data, global_entropy, global_params)
    #     self.assertEqual(res, 0.0)

if __name__ == '__main__':
    unittest.main()