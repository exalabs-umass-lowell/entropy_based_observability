import unittest
import sys

from random import Random
random = Random(1337)

sys.path.append('.')
sys.path.append('./src')

# local files
from src.setup import *

class TestSetup(unittest.TestCase):

    def test_get_initial_config(self):

        # given
        num_agents = random.randint(0, 10); num_sites = random.randint(10, 20) 
        # when
        config = initialize_sites(num_agents, num_sites)
        # then
        # should get a configuration for a given number of agents and sites
        # total elements of configuration should be equal to number of sitems
        self.assertEqual(len(config), num_sites)
        # total number of 1s should be equal to number of agents
        self.assertEqual(config.tolist().count(1), num_agents)
        # total number of -1s should be equal to number of remaining sites
        self.assertEqual(config.tolist().count(-1), num_sites - num_agents)

        # given
        num_agents = 0; num_sites = random.randint(1,100)         
        # when
        config = initialize_sites(num_agents, num_sites)
        # then
        # should contain no occupied sites when there are no agents present
        self.assertEqual((1 in config), False)



    def test_get_local_states_from_config(self):
        # given
        num_agents = random.randint(0, 10); num_sites = random.randint(10, 20)
        len_of_time = random.randint(10,60)
        # TODO this test idealy would not rely on these calls
        init_config = initialize_sites(num_agents, num_sites)
        sys_config = initialize_sites_through_time(init_config, len_of_time)
        influential_range = random.randint(1,10)
        agent_index = random.randint(0,len(init_config))
        t = random.randint(0, len_of_time)
        # when
        local_states = get_local_sites(sys_config, influential_range, agent_index, t)
        # then
        # should return an array of the states which has a length equal to the influential range
        self.assertEqual(len(local_states), influential_range)


    def test_get_hash_pose2Agent(self):
        # given
        num_agents = random.randint(0, 10); num_sites = random.randint(10, 20)
        len_of_time = random.randint(10,60)
        # TODO this test idealy would not rely on this cal
        init_config = initialize_sites(num_agents, num_sites)
        # when
        hash_pose2Agent = initialize_sites_agent_occupancy_through_time(init_config, len_of_time)
        # then
        # should 
        self.assertEqual(hash_pose2Agent.shape[0],len(init_config))
        self.assertEqual(hash_pose2Agent.shape[1],len_of_time)


    def test_get_full_config(self):
        # given
        num_agents = random.randint(0, 10); num_sites = random.randint(10, 20)
        len_of_time = random.randint(10,60)
        # TODO this test idealy would not rely on this call
        init_config = initialize_sites(num_agents, num_sites)
        # when
        full_config = initialize_sites_through_time(init_config, len_of_time)
        # then
        # should 
        self.assertEqual(full_config.shape[0],len(init_config))
        self.assertEqual(full_config.shape[1], len_of_time)

if __name__ == '__main__':
    unittest.main()