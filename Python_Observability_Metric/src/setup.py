

# dependencies - numpy, matplotlib
import numpy as np 
from numpy.random import rand
from numpy import log, dot, e
import matplotlib.pyplot as plt


# local files
from utility import *


def initialize_sites(num_of_agents, num_of_sites):
    """Create the initial sites for the traffic system
    initalized to a random configuration of [1,-1] values with
    1s representing occupied states and -1s representing unoccupied states.
    The sites are represeted as an array with a number of elements equal to 
    the given number of sites and a number of occupied
    states equal to the given number of agents

    Args:
        num_of_agents (int): the number of agents in the system
        num_of_sites (int): the number of sites in the system

    Returns:
        np.array: the array containing 1s and 0s representing occupied / unoccupied states of the traffic system
    """
    # begin with an array poplated with approprate number of ones and zeros
    initial_sites = np.append(np.ones((num_of_agents, 1)), np.zeros((num_of_sites-num_of_agents, 1)))
    # create a random permutation of the arary
    initial_sites = np.random.permutation(initial_sites)
    # scale the elemennts up by a factor of 2
    initial_sites = initial_sites*2
    # shift the elements down so they are now all 1/-1
    initial_sites = initial_sites - np.ones(num_of_sites)
    # return the resultu
    return initial_sites;


def get_local_sites(sites_through_time, influential_range, ind, t):
    """Get the local states within the influential range of a certain agent
    at a certain point in time from the system
    loops back to start position when site length is exceeded to represent
    looping road

    Args:
        sites_through_time (np.array)): 
            an array representing the system configuration over time with dimensions [num_of_sites x len_of_time]
        ind (int): 
            the index of an agent corresponding to a site in the system
        influential_range (int): 
            the range of influence for an agent ie. 
            the number of adjacent sites which will be taken account of 
            by the agent in a given site
        t (int): the current point in time

    Returns:
        np.array: an array with a number of elements equal to the influential range
            containing 1,0 values representing the occupied/unoccupied state of sites
            within that range
    """
    num_of_sites = len(sites_through_time);
    local_states = np.zeros((influential_range, 1))
    for iter in np.arange(ind, ind + influential_range):
        if(iter > num_of_sites - 1):
            #local_states[iter-ind+1] = sites_through_time[iter-num_of_sites+1, t]
            local_states[iter-ind] = sites_through_time[iter-num_of_sites, t]
        else:
            #local_states[iter-ind+1] = sites_through_time[iter, t]
            local_states[iter-ind] = sites_through_time[iter, t]
    # convert  the resulting array from 1,-1 values to 1,0 values
    local_states = convert_array_to_unity_values(local_states)
    return local_states


def initialize_sites_agent_occupancy_through_time(initial_sites, len_of_time):
    """initialize a multidimensional array which represents the positions of the agents over time
    with initial positions set
    the array will have dimension equal to [num_of_sites x len_of_time]
    

    Args:
        initial_sites (np.array):  an array representing the configuration of the initial system sites
        len_of_time (int): an integer representing the length of time of the simulation

    Returns:
        [type]: an matrix representing the positions of the agents over the duration of the simulation
            with indices corresponding weth the initial positions 
            set to values represent each agents position 
    """
    # get the number of sites in the system
    num_of_sites = len(initial_sites);
    # create an array of integers with the dimensions of system states over time
    sites_agent_occupancy_through_time = -np.ones((num_of_sites, len_of_time)).astype(int)
    # loop over the initial system configuration sites
    agent_ind = 0
    for Ind in np.arange(0, num_of_sites-1):
        # record the index of each site which corresponds with an agent ('1' for occupied site)
        if(initial_sites[Ind] == 1):
            sites_agent_occupancy_through_time[Ind, 0] = agent_ind
            agent_ind = agent_ind + 1
    return sites_agent_occupancy_through_time


def initialize_sites_through_time(initial_sites, len_of_time):
    """Get a matrix which represents the states of the system throughout the
    duration of the simulation

    Args:
        initial_sites (np.array): an array representing the initial system configuration
        len_of_time (int): an integer representing the length of time of the simulation

    Returns:
        np.array: an array of shape [num_of_sites x len_of_time] representing the configuration of the system over time
            correspinding with system states throughout the duraton of the simulation
    """
    num_of_sites = len(initial_sites);
    sites_through_time = - np.ones((num_of_sites, len_of_time))
    sites_through_time[:,0] = initial_sites
    return sites_through_time

def initialize_simulation_data(global_params):
    """ Initialize the simulation data given the global simulation parameters

    Args:
        global_params (dict): the global parameters of the simulation

    Returns:
        dict: a dictionary containing arrays which represent the simulation data over time
    """
    # extract global parameters
    num_of_sites = global_params["num_of_sites"]
    num_of_agents = global_params["num_of_agents"]
    len_of_time = global_params["len_of_time"]
    # initialize np arrays to represent the space and time data of the simulation for the sites and agents
    initial_sites = initialize_sites(num_of_agents, num_of_sites) 
    sites_through_time = initialize_sites_through_time(initial_sites, len_of_time)
    sites_agent_occupancy_through_time = initialize_sites_agent_occupancy_through_time(initial_sites, len_of_time)
    agents_motion_through_time = np.zeros((num_of_agents, len_of_time))
    # bundle site and agent data into simulation data dict
    simulation_data = dict({
        "sites":sites_through_time,
        "sites_agent_occupancy": sites_agent_occupancy_through_time,
        "agents_motion": agents_motion_through_time
    })
    # return the resulting data bundle
    return simulation_data;