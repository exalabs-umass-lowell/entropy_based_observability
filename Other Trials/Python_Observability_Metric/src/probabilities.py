# dependencies - numpy
import numpy as np 

# local files
from utility import *

# This function helps determine the probability of the local states
def calculate_probs_of_possible_local_sites(influential_range):
    """Calculate a probability distribution indicating the chances for
    each possible configuration of local sites

    Args:
        influential_range (int):
            an integer representing the number of sites over which an agent is locally influenced

    Returns:
        np.array: an array of length [2^influential_range] containing 
        equal values in the [0-1] range which represent the probabiility of possible
        local site configurations
        (this is a probability distribution so its elements will sum to 1)   
    """
    num_of_local_config = 2**influential_range
    probs_of_possible_local_sites = np.zeros((num_of_local_config, 1))
    for index in np.arange(0, num_of_local_config):
        probs_of_possible_local_sites[index] = 1/num_of_local_config
    return probs_of_possible_local_sites


# This function determines the probability at the states of being moving foward during measurements Y
def calculate_prob_ff_from_recent_motion(recent_agent_motion):
    """ Calculate the probability of free-flow motion from the
    given motion data which representing a time window view of the simulation motion history.
    The motion data consist of an array of 1 or 0 values which represent agent motion.
    The probability is calculated by applying a Gaussian kernal which weights the more recent motion values
    over the earlier ones 

    Args:
        recent_agent_motion (np.array):
            an array representing a time window slice of an agents motion durring the simulation
            the array has elements of value (1,0) and a length of [len_of_time]
            the elements represent whether a given agent transitioned to the site in front of it
            after the given moment in time

    Returns:
        float: a float between 0 and 1 representing the probability of free-flow motion from the given motion data
    """
    time_window_len = len(recent_agent_motion)
    gaussian_seq = gausswin(time_window_len*2)
    half_gaussian_seq = gaussian_seq[0:time_window_len]
    # gaussianKernel = half_gaussian_seq / (np.sum(gaussian_seq)*2)
    gaussianKernel = half_gaussian_seq / np.sum(half_gaussian_seq)
    prob_ff_from_recent_motion = np.dot(gaussianKernel.transpose(), recent_agent_motion)
    return prob_ff_from_recent_motion

 
 
def calculate_prob_ff_from_local_site_config(local_sites, global_params):
    """Calculate the probability for a free-flow state given a specific configuration
    of local sites

    Args:
        local_sites (np.array): an array with a number of length [influential_range]
            containing 1,0 values representing the occupied/unoccupied state of sites
            within that range
        global_params (dict): a dictionary containing the global simulation parameters

    Returns:
        float: a floating point number representing the probability for free-flow
    """
    external_field_coeff = global_params["external_field_coeff"]
    weighted_local_interactions = calculate_weighted_local_interactions(global_params)
    inter_coeff_con_trans = weighted_local_interactions.conjugate().transpose()
    coeff_prod = np.dot(inter_coeff_con_trans, local_sites.flatten()) 
    expFact = np.exp(external_field_coeff - coeff_prod)
    prob = min(1.0, expFact);
    return prob
        

def calculate_prob_dist_agents_in_sites(num_of_local_agents, influential_range):
    """Calculate a probability distribution indicating the chances for
    each possible way in which the given number of agents can be fit into the sites
    which fall in the influential range

    Args:
        num_of_local_agents (int): the number of agents which are found within the influential range
        influential_range (int): an integer representing the number of sites over which an agent is influenced

    Returns:
        np.array: an array of length [nchoosek(influential_rangeo num_of_local_agents)] containing
            equal values in the [0-1] range which represent the probability of a choice of agent
            placement within the influential range. 
    """
    num_of_possible_local_sites = nchoosek(influential_range, num_of_local_agents);
    prob_of_possible_local_sites = 1 / num_of_possible_local_sites;
    res = prob_of_possible_local_sites * np.ones(num_of_possible_local_sites).transpose()
    return res


# TODO FOUND MATCHING original MATLAB
def calculate_prob_dist_agents_out_of_range(num_of_agents_out_of_range, num_of_sites_out_of_range):
    """Calculate a probability distribution for the agents that 
    fall outside of the influential range
    TODO ????????????????????? need better understanding of this

    Args:
        num_of_agents_out_of_range (int): the number of agents which are found in sites outside of the influential range
        num_of_sites_out_of_range (int): the number of sites which are outside of the influential range

    Returns:
        np.array: an array of length [num_of_agents - 1 - num_local_agents]
            representing a probability distribution for the out of range agents
            (being a distrubution, it contains probabiliy values that sum to 1)
    """
    energy_dist_other_agents = calculate_agents_energy_level_count(num_of_agents_out_of_range, num_of_sites_out_of_range);
    res = calculate_prob_dist(energy_dist_other_agents)
    return res

def calculate_energy_distribution_other_sites(num_of_agents, num_of_sites):
    """[summary]

    Args:
        num_of_agents (int): the number of agents in the system
        num_of_sites (int): the number of sites in the system

    Returns:
        [type]: [description]
    """
    # res = (num_of_sites - 2*num_of_agents - 1) -  2*np.arange(num_of_agents, num_of_sites-1, 4) - 2;
    start_val = (num_of_sites - 2*num_of_agents - 1) - 2*num_of_agents
    end_val = num_of_sites-1 - 2
    step = 4
    res = np.arange(start_val, end_val, step);
    return res


# This function helps determine the probability of the global states
def calculate_agents_energy_level_count(num_of_agents, num_of_sites):
    """[summary]
    TODO ?????????????????? need better understanding of this

    Args:
        num_of_agents (int): the number of agents in the system
        num_of_sites (int): the number of sites in the system

    Returns:
        np.array: an array of length [num_of_agents] containing values which represent ???
    """
    # TODO figure out how this energy level is derived
    energy_distrubution_other_sites = calculate_energy_distribution_other_sites(num_of_agents, num_of_sites)
    num_of_level = energy_distrubution_other_sites.shape[0]
    num_of_vacant_sites = num_of_sites - num_of_agents
    count_of_energy_level = np.zeros(num_of_level)
    # count_of_energy_level = np.zeros(num_of_agents)
    # find spaces between vacant sites and insert set of vehicles  *  group vehicles into (index-1) sets
    # for index in np.arange(0, num_of_level):
    for index in np.arange(0, num_of_agents):
        num_of_counted_agents = index + 1
        num_of_remaining_agents = num_of_agents-num_of_counted_agents
        # if(num_of_vacant_sites >= num_of_agents+1-index):
        if(num_of_vacant_sites > num_of_remaining_agents):
            num_of_ways_to_put_remaining_agents_in_vacant_sites = nchoosek(num_of_vacant_sites, num_of_remaining_agents)
            num_of_ways_to_select_counted_agents = nchoosek(num_of_agents, num_of_counted_agents)
            # count_of_energy_level[index-1] = nchoosek(num_of_vacant_sites-1, num_of_agents+1-index)*nchoosek(num_of_agents-1, index-1);
            count_of_energy_level[index] = num_of_ways_to_put_remaining_agents_in_vacant_sites * num_of_ways_to_select_counted_agents
        else:
            num_of_ways_to_fill_vacent_sites_with_remaining_agents = nchoosek(num_of_remaining_agents, num_of_vacant_sites)
            if(num_of_counted_agents < num_of_vacant_sites):
                num_of_ways_to_have_selected_vacent_sites_with_counted_agents = nchoosek(num_of_vacant_sites, num_of_counted_agents)
            else:
                num_of_ways_to_have_selected_vacent_sites_with_counted_agents = 0
            # count_of_energy_level[index-1] = nchoosek(num_of_agents+1-index, num_of_vacant_sites-1)*nchoosek(num_of_vacant_sites-1, index-1);
            count_of_energy_level[index] = num_of_ways_to_fill_vacent_sites_with_remaining_agents * num_of_ways_to_have_selected_vacent_sites_with_counted_agents;
    return count_of_energy_level

    
def calculate_prob_dict_local_sites_given_motion_and_agent_num(prob_dict_local_sites_given_motion, num_agents, global_params):
    possible_local_sites = calculate_possible_sites(global_params["influential_range"])
    num_remote_sites = global_params["num_of_sites"] - global_params["influential_range"]
    min_local_agents = max(num_agents - num_remote_sites - 1, 0)
    max_local_agents = min(num_agents, global_params["influential_range"])
    # reduce the prob_dict to only sites possible with the given number of agents
    prob_dict_local_sites_given_motion_and_agent_num = prob_dict_local_sites_given_motion.copy()
    for local_sites in possible_local_sites:
        local_agent_count = list(local_sites).count(1)
        if local_agent_count > max_local_agents or local_agent_count < min_local_agents:
            prob_dict_local_sites_given_motion_and_agent_num[str(local_sites)] = 0.0
    norm_sum = np.sum(list(prob_dict_local_sites_given_motion_and_agent_num.values()))
    for key in prob_dict_local_sites_given_motion_and_agent_num.keys():
        val = prob_dict_local_sites_given_motion_and_agent_num[key]
        prob_dict_local_sites_given_motion_and_agent_num[key] = val / norm_sum
    return prob_dict_local_sites_given_motion_and_agent_num


def calculate_prob_of_agent_num(num_of_agents, num_of_sites):
    num_of_non_agent_sites = num_of_sites - 1
    possible_sites = calculate_possible_sites(num_of_non_agent_sites)
    sites_with_agent_num = 0
    for sites in possible_sites:
        if list(sites).count(1) == num_of_agents - 1:
            sites_with_agent_num = sites_with_agent_num + 1
    prob_of_agent_num = sites_with_agent_num / len(possible_sites)
    return prob_of_agent_num


def calculate_prob_agent_num_given_local_sites(num_of_agents, num_of_sites, local_sites):
    num_of_local_sites = len(local_sites)
    num_of_remote_sites = num_of_sites - num_of_local_sites - 1
    num_of_local_agents = list(local_sites).count(1)
    num_of_remote_agents = num_of_agents - num_of_local_agents - 1
    possible_remote_sites = calculate_possible_sites(num_of_remote_sites)
    if 1 + num_of_local_agents + num_of_remote_sites < num_of_agents:
        return 0.0
    elif num_of_remote_agents < 0:
        return 0.0
    else:
        num_of_possible_remote_sites_given_agent_num = nchoosek(num_of_remote_sites, num_of_remote_agents)
        num_of_possible_remote_sites = len(possible_remote_sites)
        return num_of_possible_remote_sites_given_agent_num / num_of_possible_remote_sites 

        

def calculate_prob_of_local_agent_num(num_local_agents, possible_local_sites):
    site_count = 0
    total_sites = len(possible_local_sites)
    for local_sites in possible_local_sites:
        if list(local_sites).count(1) == num_local_agents:
            site_count = site_count + 1
    return site_count / total_sites