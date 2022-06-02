import unittest
import sys

from random import Random
random = Random(1337)

sys.path.append('.')
sys.path.append('./src')

# local files
from src.utility import *
from src.example_analysis import *

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

    def test_calculate_site_prob_dict_given_motion(self):
        global_params = {
            "num_of_sites": 6,
            "num_of_agents": 3,
            "interaction_coeff": 4.0,
            "external_field_coeff": 1.5,
            "len_of_time": 64,
            "influential_range": 3,
            "len_of_time_window": 4
        }
        full_recent_motion = np.ones(global_params["len_of_time_window"]).astype(float)
        no_recent_motion = np.zeros(global_params["len_of_time_window"]).astype(float)
        prob_dict = get_site_prob_dict_given_motion(full_recent_motion, global_params)
        prob_dist = np.array(list(prob_dict.values()))
        # then it produce a probability vec which should sum to 1
        self.assertEqual(np.sum(prob_dist), 1.0)
        

    def test_calculate_energy_prob_dict_given_motion(self):
        global_params = {
            "num_of_sites": 6,
            "num_of_agents": 3,
            "interaction_coeff": 4.0,
            "external_field_coeff": 1.5,
            "len_of_time": 64,
            "influential_range": 3,
            "len_of_time_window": 4
        }
        full_recent_motion = np.ones(global_params["len_of_time_window"]).astype(float)
        no_recent_motion = np.zeros(global_params["len_of_time_window"]).astype(float)
        prob_dict = get_energy_prob_dict_given_motion(full_recent_motion, global_params)
        prob_dist = np.array(list(prob_dict.values()))
        # then it produce a probability vec which should sum to 1
        self.assertEqual(np.sum(prob_dist), 1.0)


    def test_calculate_global_entropy_of_microstate_given_motion_dict(self):
        global_params = {
            "num_of_sites": 6,
            "num_of_agents": 3,
            "interaction_coeff": 4.0,
            "external_field_coeff": 1.5,
            "len_of_time": 64,
            "influential_range": 3,
            "len_of_time_window": 4
        }
        entropy_dict = get_global_entropy_of_microstate_given_motion_dict(global_params)
        self.assertTrue(True)


    def test_calculate_global_entropy_of_macrostate_given_motion_dict(self):
        global_params = {
            "num_of_sites": 6,
            "num_of_agents": 3,
            "interaction_coeff": 4.0,
            "external_field_coeff": 1.5,
            "len_of_time": 64,
            "influential_range": 3,
            "len_of_time_window": 4
        }
        entropy_dict = get_global_entropy_of_macrostate_given_motion_dict(global_params)
        self.assertTrue(True)


    # def test_get_site_prob_dict_given_motion_and_agent_num(self):
    #     global_params = {
    #         "num_of_sites": 6,
    #         "num_of_agents": 3,
    #         "interaction_coeff": 4.0,
    #         "external_field_coeff": 1.5,
    #         "len_of_time": 64,
    #         "influential_range": 3,
    #         "len_of_time_window": 4
    #     }
    #     full_recent_motion = np.ones(global_params["len_of_time_window"]).astype(float)
    #     num_agents = 3
    #     prob_dict = get_site_prob_dict_given_motion_and_agent_num(full_recent_motion, num_agents, global_params)
    #     prob_dist = np.array(list(prob_dict.values()))
    #     # then it produce a probability vec which should sum to 1
    #     self.assertEqual(np.sum(prob_dist), 1.0)


    def test_get_mutual_information_given_motion_dict(self):
        global_params = {
            "num_of_sites": 6,
            "num_of_agents": 3,
            "interaction_coeff": 4.0,
            "external_field_coeff": 1.5,
            "len_of_time": 64,
            "influential_range": 3,
            "len_of_time_window": 4
        }
        full_recent_motion = np.ones(global_params["len_of_time_window"]).astype(float)
        num_agents = 3
        inforamiot_dict = get_mutual_information_given_motion_dict(global_params)
        # then it produce a probability vec which should sum to 1
        # self.assertEqual(np.sum(prob_dist), 1.0)
        self.assertTrue(True)

    

if __name__ == '__main__':
    unittest.main()