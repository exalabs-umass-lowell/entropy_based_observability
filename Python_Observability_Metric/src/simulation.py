# dependencies - numpy, matplotlib
import numpy as np 
from numpy.random import rand
from numpy import log, dot, e
import matplotlib.pyplot as plt


from utility import *
from setup import *

def calculate_neighboring_indices_given_boundary(i, num_of_sites):
    """ Calculate the pair of site indices which should be considered for possible
    transitions given the index of a site and teh total number of sites
    (this accounts for the boundary of the simulation space by looping the n+1th site back to the 0th)

    Args:
        i (int): index of a given site
        num_of_sites (int): total number of sites in the system

    Returns:
        int, int: integers representing the indices of the two sites to be considered
    """
    ind = i
    ind2 = i+1
    if(ind == num_of_sites - 1):
        ind2 = 0
    return ind, ind2

def calculate_transmission_prob(sites_through_time, global_params, ind, t):
    """Calculate the probability for a given site to transition with the site 
    in front of it (ie: simulate one step of a car moving through traffic)

    Args:
        sites_through_time (np.array): 
            an array with shape total_site_num x len_of_time which represents the site states throughout the simulation
        global_params (dict): 
            a dictionary containing the global parameters of the simulation
        ind (int)): 
            an integer representing the index of the site in question
        t (int): 
            an integer representing the time in question

    Returns:
        float: the probability of the site to transition with the site in front of it
    """
    H = calculate_effective_hamiltonion(sites_through_time, global_params, ind, t)
    p_tr = get_prob_from_hamiltonion(H)
    return p_tr

def is_valid_transition(sites_through_time, ind, ind2, t):
    """Check  if a given pair of sites are valid for a site transition to occur

    Args:
        sites_through_time (np.array): 
            an array with shape total_site_num x len_of_time which represents the site states throughout the simulation
        ind (int)): 
            an integer representing the index of one of the sites in question
        ind2 (int):
            an integer representing the index of the other site in question
        t (int): 
            an integer representing the time in question

    Returns:
        bool: a boolean value representing whether or not the first site can transition to the next site
    """
    # get the states corresponding with the sites in question
    site1 = sites_through_time[ind, t]
    site2 = sites_through_time[ind2, t]
    # if the first site is occupied and the other is empty then a transition can occur
    result = site1 == 1 and site2 == -1
    return result

def record_successful_transition(simulation_data, ind, ind2, t):
    """Record the transition of an agent at a given time in the simulation
    by changing the states of sites in the simulation data

    Args:
        simulation_data (dict): 
            a dictionary containg the np.arays which keep track of the simulation site and agent data
        ind (int)): 
            an integer representing the index of one of the sites in question
        ind2 (int):
            an integer representing the index of the other site in question
        t (int): 
            an integer representing the time in question

    Returns:
        dict: the simulation data dictionary now with arrays containing modified states
    """
    # extract simulation data
    sites_through_time = simulation_data["sites"]
    sites_agent_occupancy_through_time = simulation_data["sites_agent_occupancy"]
    agents_motion_through_time = simulation_data["agents_motion"]
    # swap 
    sites_through_time[ind, t + 1] = sites_through_time[ind2, t]
    sites_through_time[ind2, t + 1] = sites_through_time[ind, t]
    # record the agent motion in the agent data 
    agent_ind = sites_agent_occupancy_through_time[ind, t]
    agents_motion_through_time[agent_ind, t + 1] = 1
    sites_agent_occupancy_through_time[ind2, t + 1] = agent_ind
    # save the values back into the bundled simulation data
    simulation_data["sites"] = sites_through_time
    simulation_data["sites_agent_occupancy"] = sites_agent_occupancy_through_time 
    simulation_data["agents_motion"] = agents_motion_through_time 
    return simulation_data

def record_failed_transition(simulation_data, ind, ind2, t):
    """Record the failure of a transition of an agent at a given time in the simulation
    by saving the states of sites in the simulation data

    Args:
        simulation_data (dict): 
            a dictionary containg the np.arays which keep track of the simulation site and agent data
        ind (int)): 
            an integer representing the index of one of the sites in question
        ind2 (int):
            an integer representing the index of the other site in question
        t (int): 
            an integer representing the time in question

    Returns:
        dict: the simulation data dictionary now with arrays containing modified states
    """
    # extract simulation data
    sites_through_time = simulation_data["sites"]
    sites_agent_occupancy_through_time = simulation_data["sites_agent_occupancy"]
    # save the same site value and agent index value for the site and site agent occupancy at the next time step
    sites_through_time[ind, t + 1] = sites_through_time[ind, t]
    sites_agent_occupancy_through_time[ind, t + 1] = sites_agent_occupancy_through_time[ind, t]
    # save the values back into the bundled simulation data
    simulation_data["sites"] = sites_through_time
    simulation_data["sites_agent_occupancy"] = sites_agent_occupancy_through_time 
    return simulation_data

def attempt_site_transition(simulation_data, global_params, ind, ind2, t):
    """Attempt to perform an agent transition for a given pair of sites at a given time in the simulation
    record the resulting states of these sites given either a success or failure for that time in the simulation data

    Args:
        simulation_data (dict): 
            a dictionary containg the np.arays which keep track of the simulation site and agent data
        global_params (dict): 
            a dictionary containing the global parameters of the simulation
        ind (int)): 
            an integer representing the index of one of the sites in question
        ind2 (int):
            an integer representing the index of the other site in question
        t (int): 
            an integer representing the time in question

    Returns:
        dict, int: 
            the simulation data dictionary now with arrays containing modified states, and the possiblymodified site index 
            the site index must be returned since a successful site transation will increment it by one
    """
    # extract simulation data
    sites_through_time = simulation_data["sites"]
    # % Metroplis Algorithm to update the sites
    # calculate the transmission probability and a random value to see if a transmission has occured
    P_tr = calculate_transmission_prob(sites_through_time, global_params, ind, t)
    rand_val = np.random.rand()             # produce a random number between 0 and 1
    # make sure the vehicle are moving forward/backward
    if(is_valid_transition(sites_through_time, ind, ind2, t) and rand_val <= P_tr):
        simulation_data = record_successful_transition(simulation_data, ind, ind2, t)
        ind = ind + 1
    else:
        simulation_data = record_failed_transition(simulation_data, ind, ind2, t)
    # return the modified data and site index
    return simulation_data, ind

def step_configuration(simulation_data, global_params, t):
    """Perform the operations necessary to simulate one time step of the simulation
    save the resulting changes to the site configuration to the simulation data

    Args:
        simulation_data (dict): 
            a dictionary containg the np.arays which keep track of the simulation site and agent data
        global_params (dict): 
            a dictionary containing the global parameters of the simulation
        t (int): 
            an integer representing the time in question

    Returns:
        dict: the simulation data dictionary now with arrays containing modified states
    """
    num_of_sites = global_params["num_of_sites"]
    i = 0
    while(i <= num_of_sites - 1):
        num_of_sites = global_params["num_of_sites"]
        # In case get to the boundary
        ind, ind2 = calculate_neighboring_indices_given_boundary(i, num_of_sites);
        simulation_data, i = attempt_site_transition(simulation_data, global_params, ind, ind2, t);
        # increment the site counter
        i = i + 1
        # return the simulation data
    return simulation_data

def get_agent_start_positions(sites_agent_occupancy_through_time, global_params):
    """ Get an array which represents the starting positions for each agent in the simulation

    Args:
        sites_agent_occupancy_through_time (np.array): 
            thearray representing the sites during the simulation
            with dimensions the num_of_site x len_of_tme
        global_params (dict): 
            global parameters used throughout the simulation

    Returns:
        np.array: 
            an array containing initial site indices for each agent
    """
    num_of_agents = global_params["num_of_agents"]
    agent_start_positions = np.zeros((num_of_agents, 1))
    # for i in range(0, sites_agent_occupancy_through_time.shape[1]):
    site = sites_agent_occupancy_through_time[:,0]
    for i in range(0, len(site)):
        occupant = site[i]
        if(occupant != -1):
            agent_start_positions[occupant] = i
    return agent_start_positions

def calculate_effective_hamiltonion(sites_through_time, global_params, ind, t):
    """Calculate the total ammount of energy for a given time and place,
    ie: the effective Hamiltonion given the states of local sites
    around a certain site at a certain time

    Args:
        sites_through_time (np.array): 
            an array with shape total_site_num x len_of_time which represents the site states throughout the simulation
        global_params (dict): 
            a dictionary containing the global parameters of the simulation
        ind (int)): 
            an integer representing the index of the site in question
        t (int): 
            an integer representing the time in question

    Returns:
        float: a floating point value representing the total energy of the system at the given time and place
            ie: the effective Hamiltonion of the system
    """
    influential_range = global_params["influential_range"]
    # Interaction terms
    site = sites_through_time[ind, t]
    local_sites = get_local_sites(sites_through_time, influential_range, ind, t)
    H = get_hamiltonion_for_site_from_local_sites(site, local_sites, global_params)
    return H

def perform_simulation(simulation_data, global_params):
    """Perform the tasks to simulate each step of the simulation
    for the full duration of time

    Args:
        simulation_data (dict): 
            a dictionary containg the np.arays which keep track of the simulation site and agent data
        global_params (dict): 
            a dictionary containing the global parameters of the simulation

    Returns:
        dict: the simulation data dictionary now with arrays containing 
            the site and agent data of the simulation
    """
    # extract global params
    len_of_time = global_params["len_of_time"]
    # perform a step for each tick of the duration of time in the simulation
    for t in np.arange(0, len_of_time - 1):
        simulation_data = step_configuration(simulation_data, global_params, t)
    return simulation_data


def statistical_mechanics_based_traffic(global_params):
    """ Setup the initial site configuration with agents, then run the simulation
    recording the results in the simulation data which is finally returned

    Args:
        global_params (dict): 
            global parameters used throughout the simulation

    Returns:
        dict: the simulation data dictionary now with arrays containing 
            the site and agent data of the simulation
    """
    # initialize site and agent simulation data
    simulation_data = initialize_simulation_data(global_params);
    # run the simulation function on the data
    simulation_data = perform_simulation(simulation_data, global_params)
    # return the result
    return simulation_data
