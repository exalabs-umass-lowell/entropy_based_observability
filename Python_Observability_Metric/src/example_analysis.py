# dependencies - numpy
from glob import glob
import numpy as np 

# local files
from probabilities import *
from analysis import *
from utility import *


def get_site_prob_given_motion_dict___with_know_agent_num(motion, global_params):
    possible_local_sites = calculate_possible_sites(global_params["influential_range"])
    prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(motion, global_params)
    # prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(no_recent_motion, global_params)
    x_data = [str(local_sites) for local_sites in possible_local_sites]
    y_data = [item.item() for item in prob_dist_local_sites_given_motion]
    prob_dict_local_sites_given_motion = {}
    for A, B in zip(x_data, y_data):
        prob_dict_local_sites_given_motion[str(A)] = B

    prob_sum = np.sum(prob_dist_local_sites_given_motion)
    if  (1.0 - prob_sum) > 0.001:
        print("PROBLEM")

    site_prob_given_motion_dict = dict()
    num_agents = global_params["num_of_agents"]
    # max_agents = global_params["num_of_sites"]
    num_of_remote_sites = global_params["num_of_sites"] - global_params["influential_range"] - 1
    # possible_num_of_agent = list(range(0, max_agents + 1))
    # for num_agents in possible_num_of_agent:
    for local_sites in possible_local_sites: 
        prob_local_sites_given_motion = prob_dict_local_sites_given_motion[str(local_sites)]
        possible_remote_sites, prob_dist = calculate_prob_dist_remote_sites_given_local_sites_and_agent_num(local_sites, num_agents, global_params)
        # make prob_dist uniform since agent number is unknown
        # prob_dist[prob_dist >= 0] = 1 / len(possible_remote_sites)
        prob_dict_remote_sites_given_local = {}
        for A, B in zip(possible_remote_sites, prob_dist):
            prob_dict_remote_sites_given_local[str(A)] = B
        if np.sum(prob_dist) != 1.0:
            continue
        for remote_sites in possible_remote_sites:
            constructed_sites = np.concatenate((local_sites, remote_sites), axis=None)
            # if sites_invalid_for_agent_number(local_sites, remote_sites, num_agents):
            #     # site_prob_given_motion_dict[str(constructed_sites)] = 0.0
            #     continue
            # else:
            prob_remote_sites_given_local = prob_dict_remote_sites_given_local[str(remote_sites)]
            prob_global_sites_given_motion = prob_remote_sites_given_local * prob_local_sites_given_motion
            if str(constructed_sites) not in site_prob_given_motion_dict:
                site_prob_given_motion_dict[str(constructed_sites)] = prob_global_sites_given_motion
            else: 
                site_prob_given_motion_dict[str(constructed_sites)] = site_prob_given_motion_dict[str(constructed_sites)] + prob_global_sites_given_motion
    prob_sum = np.sum(list(site_prob_given_motion_dict.values()))
    if  (1.0 - prob_sum) > 0.001 :
        print("PROBLEM site_prob_given_motion comes to: " + sum(prob_sum))
    return site_prob_given_motion_dict


def get_site_prob_given_motion_dict(motion, global_params):
    possible_local_sites = calculate_possible_sites(global_params["influential_range"])
    prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(motion, global_params)
    # prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(no_recent_motion, global_params)
    x_data = [str(local_sites) for local_sites in possible_local_sites]
    y_data = [item.item() for item in prob_dist_local_sites_given_motion]
    prob_dict_local_sites_given_motion = {}
    for A, B in zip(x_data, y_data):
        prob_dict_local_sites_given_motion[str(A)] = B

    prob_sum = np.sum(prob_dist_local_sites_given_motion)
    if  (1.0 - prob_sum) > 0.001:
        print("PROBLEM")

    site_prob_given_motion_dict = dict()
    # num_agents = global_params["num_of_agents"]
    max_agents = global_params["num_of_sites"]
    num_of_remote_sites = global_params["num_of_sites"] - global_params["influential_range"] - 1
    # possible_num_of_agent = list(range(0, max_agents + 1))
    # for num_agents in possible_num_of_agent:
    for local_sites in possible_local_sites: 
        prob_local_sites_given_motion = prob_dict_local_sites_given_motion[str(local_sites)]
        possible_remote_sites, prob_dist = calculate_prob_dist_remote_sites_given_local_sites_and_agent_num(local_sites, max_agents, global_params)
        # make prob_dist uniform since agent number is unknown
        prob_dist[prob_dist >= 0] = 1 / len(possible_remote_sites)
        prob_dict_remote_sites_given_local = {}
        for A, B in zip(possible_remote_sites, prob_dist):
            prob_dict_remote_sites_given_local[str(A)] = B
        if np.sum(prob_dist) != 1.0:
            continue
        for remote_sites in possible_remote_sites:
            constructed_sites = np.concatenate((local_sites, remote_sites), axis=None)
            # if sites_invalid_for_agent_number(local_sites, remote_sites, num_agents):
            #     # site_prob_given_motion_dict[str(constructed_sites)] = 0.0
            #     continue
            # else:
            prob_remote_sites_given_local = prob_dict_remote_sites_given_local[str(remote_sites)]
            prob_global_sites_given_motion = prob_remote_sites_given_local * prob_local_sites_given_motion
            if str(constructed_sites) not in site_prob_given_motion_dict:
                site_prob_given_motion_dict[str(constructed_sites)] = prob_global_sites_given_motion
            else: 
                site_prob_given_motion_dict[str(constructed_sites)] = site_prob_given_motion_dict[str(constructed_sites)] + prob_global_sites_given_motion
    prob_sum = np.sum(list(site_prob_given_motion_dict.values()))
    if  (1.0 - prob_sum) > 0.001 :
        print("PROBLEM site_prob_given_motion comes to: " + sum(prob_sum))
    return site_prob_given_motion_dict
    
def get_array_from_sites_str(sites_str):
    result = np.array([])
    for char in sites_str.strip('[] ').replace(',',' ').split(' '):
        result = np.append(result, float(char))
    return result
        
    
def get_energy_prob_dict_given_motion(motion, global_params):
    energy_prob_dict = dict()
    site_prob_dict = get_site_prob_given_motion_dict(motion, global_params)
    for sites in site_prob_dict.keys():
        # sites_array = np.array(str(sites)) 
        sites_array = get_array_from_sites_str(sites)
        sites_prob = site_prob_dict[sites]
        sites_energy = round(get_hamiltonion_for_system_sites(sites_array, global_params), 2)
        if str(sites_energy) not in energy_prob_dict:
            energy_prob_dict[str(sites_energy)] = sites_prob
        else: 
            energy_prob_dict[str(sites_energy)] = energy_prob_dict[str(sites_energy)] + sites_prob
    return energy_prob_dict


def get_possible_sites_energy_dict(global_params):
    energy_dict = dict()
    possible_sites = calculate_possible_sites(global_params["num_of_sites"])
    for sites in possible_sites:
        sites_energy = round(get_hamiltonion_for_system_sites(sites, global_params), 2)
        if str(sites) not in energy_dict:
            energy_dict[str(sites)] = sites_energy
    return energy_dict


def get_global_entropy_of_microstate_given_motion_dict(global_params):
    possible_motion = calculate_possible_sites(global_params["len_of_time_window"])
    global_entropy_given_motion = dict()

    for motion in possible_motion:
        possible_site_probs = get_site_prob_given_motion_dict(motion, global_params)
        prob_dist = np.array(list(possible_site_probs.values()))

        entropy = calculate_entropy_from_prob_dist(prob_dist)
        global_entropy_given_motion[str(motion)] = entropy
    return global_entropy_given_motion
    
def get_global_entropy_of_macrostate_given_motion_dict(global_params):
    possible_motion = calculate_possible_sites(global_params["len_of_time_window"])
    global_entropy_given_motion = dict()
    # weighted_global_entropy_given_motion = dict()

    for motion in possible_motion:
        # possible_site_probs = get_site_prob_given_motion_dict(motion, global_params)
        possible_energy_probs = get_energy_prob_dict_given_motion(motion, global_params)
        prob_dist = np.array(list(possible_energy_probs.values()))

        entropy = calculate_entropy_from_prob_dist(prob_dist)
        global_entropy_given_motion[str(motion)] = entropy
        # also calculate this entropy weighted by the prob of ff for the next graph
        prob_motion = calculate_prob_ff_from_recent_motion(motion)
        # weighted_global_entropy_given_motion[str(motion)] = entropy * prob_motion
    return global_entropy_given_motion

    
def get_global_entropy_of_macrostate_diff_given_motion_dict(global_params):
    possible_motion = calculate_possible_sites(global_params["len_of_time_window"])
    cropped_global_entropy_given_motion = dict()
    # weighted_global_entropy_given_motion = dict()
    for motion in possible_motion:
        # possible_site_probs = get_site_prob_given_motion_dict(motion, global_params)
        possible_energy_probs = get_energy_prob_dict_given_motion(motion, global_params)
        prob_dist = np.array(list(possible_energy_probs.values()))

        entropy = calculate_entropy_from_prob_dist(prob_dist)
        cropped_global_entropy_given_motion[str(motion)] = entropy
    cropped_global_entropy_given_motion = crop_dict_values(cropped_global_entropy_given_motion)
    return cropped_global_entropy_given_motion

def crop_dict_values(dict):
    min_dict_value = min(dict.values())
    cropped_dict = {k: v - min_dict_value for k, v in dict.items()}
    return cropped_dict


def get_total_entropy_of_system_macrostate(global_params):
    energy_dict = get_possible_sites_energy_dict(global_params)
    energy_prob_dist = calculate_prob_dist(list(energy_dict.values()))
    total_entropy = calculate_entropy_from_prob_dist(energy_prob_dist)
    return total_entropy

def sites_invalid_for_agent_number(local_sites, remote_sites, agent_num):
    num_local_agents = list(local_sites).count(1)
    num_remote_agents = list(remote_sites).count(1)
    is_invalid = num_local_agents + num_remote_agents + 1 != agent_num
    return is_invalid


# def calculate_likelihood_local_sites_given_agent_motion_and_agent_number(recent_agent_motion, agent_num, global_params):
#     influential_range = global_params["influential_range"]
#     # calculate the probabilities of all possible local site configurations for the influential range
#     possible_local_sites = calculate_possible_sites(influential_range)
    
#     feasable_local_sites = np.array([])
#     for sites in possible_local_sites:
#         num_local_agents = list(sites).count(1)
#         if(num_local_agents + 1 < agent_num):
#             feasable_local_sites = np.append(feasable_local_sites, sites) 
#     probs_of_possible_local_sites  = calculate_probs_of_possible_local_sites(influential_range)
#     # calculate the probability of free flow from this
#     prob_ff_from_recent_motion = calculate_prob_ff_from_recent_motion(recent_agent_motion.transpose())
#     prob_jam_from_recent_motion = 1 - prob_ff_from_recent_motion
#     # calculate the probabilitys for possible local sites given the motion data from the time window
#     probs_ff_given_possible_local_sites = calculate_probs_ff_given_possible_local_sites(feasable_local_sites, global_params);
#     probs_jam_given_possible_local_sites = (1 - probs_ff_given_possible_local_sites)
#     # calculate the probabilities for jam / free-flow given any of the possible local sites
#     probs_ff_and_local_sites = probs_ff_given_possible_local_sites * probs_of_possible_local_sites
#     probs_jam_and_local_sites = probs_jam_given_possible_local_sites * probs_of_possible_local_sites
#     # from thes probabilities calculate the likelyhood for the local sites to be in any configuration
#     probs_ff_and_local_sites_given_agent_motion = 0
#     probs_jam_and_local_sites_given_agent_motion = 0
#     # avoid dividing by 0 if any of the probabilities comes to that
#     if(prob_ff_from_recent_motion != 0):
#         probs_ff_and_local_sites_given_agent_motion = probs_ff_and_local_sites / prob_ff_from_recent_motion
#     if(prob_jam_from_recent_motion != 0):
#         probs_jam_and_local_sites_given_agent_motion = probs_jam_and_local_sites / prob_jam_from_recent_motion
#     # the overall probability of the local site configuration is the combination of probabilities for free-flow and for jam
#     likelihood_local_sites_given_agent_motion = probs_ff_and_local_sites_given_agent_motion + probs_jam_and_local_sites_given_agent_motion
#     return likelihood_local_sites_given_agent_motion


# # derive_and_normalize_likelihood_from_bayesian
# def calculate_prob_dist_local_sites_given_motion_and_agent_num(recent_agent_motion, agent_num, global_params):
#     likelihood_local_sites_given_agent_motion = calculate_likelihood_local_sites_given_agent_motion_and_agent_number(recent_agent_motion, agent_num, global_params)
#     # normalise the likelihood arary to get a probability distribution for the result
#     probs_local_sites_given_ff = calculate_prob_dist(likelihood_local_sites_given_agent_motion)
#     return probs_local_sites_given_ff


    
def get_site_prob_dict_given_motion_and_agent_num(motion, num_agents, global_params):
    possible_local_sites = calculate_possible_sites(global_params["influential_range"])
    prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(motion, global_params)
    sum_of_prob_dist = np.sum(prob_dist_local_sites_given_motion)
    x_data = [str(local_sites) for local_sites in possible_local_sites]
    y_data = [item.item() for item in prob_dist_local_sites_given_motion]
    prob_dict_local_sites_given_motion = {}
    for A, B in zip(x_data, y_data):
        prob_dict_local_sites_given_motion[str(A)] = B

    # constructed_possible_sites = dict()
    # for local_sites in possible_local_sites: 
    #     prob_local_sites_given_motion = prob_dict_local_sites_given_motion[str(local_sites)]
    #     possible_remote_sites, prob_dist = calculate_prob_dist_remote_sites_given_local_sites_and_agent_num(local_sites, num_agents, global_params)
    #     # make prob_dist uniform since agent number is unknown
    #     # prob_dist[prob_dist >= 0] = 0.25
    #     prob_dict_remote_sites_given_local = {}
    #     for A, B in zip(possible_remote_sites, prob_dist):
    #         prob_dict_remote_sites_given_local[str(A)] = B
    #     # if np.sum(prob_dist) != 1.0:
    #         # continue
    #     for remote_sites in possible_remote_sites:
    #         constructed_sites = np.concatenate((local_sites, remote_sites), axis=None)
    #         if sites_invalid_for_agent_number(local_sites, remote_sites, num_agents):
    #             # continue
    #             constructed_possible_sites[str(constructed_sites)] = 0.0
    #         else:
    #             prob_remote_sites_given_local = prob_dict_remote_sites_given_local[str(remote_sites)]
    #             prob_global_sites_given_motion = prob_remote_sites_given_local * prob_local_sites_given_motion
    #             if str(constructed_sites) not in constructed_possible_sites:
    #                 constructed_possible_sites[str(constructed_sites)] = prob_global_sites_given_motion
    #             else: 
    #                 constructed_possible_sites[str(constructed_sites)] = constructed_possible_sites[str(constructed_sites)] + prob_global_sites_given_motion
    # return constructed_possible_sites


    constructed_possible_sites = dict()
    # for num_agents in range (0,6):
    prob_dict_local_sites_given_motion_and_agent_num = calculate_prob_dict_local_sites_given_motion_and_agent_num(prob_dict_local_sites_given_motion, num_agents, global_params)
    final_sum = np.sum(list(prob_dict_local_sites_given_motion_and_agent_num.values()))
    for local_sites in possible_local_sites: 
        prob_local_sites_given_motion_and_agent_num = prob_dict_local_sites_given_motion_and_agent_num[str(local_sites)]
        possible_remote_sites, prob_dist = calculate_prob_dist_remote_sites_given_local_sites_and_agent_num(local_sites, num_agents, global_params)
        sum = np.sum(prob_dist)
        prob_dict_remote_sites_given_local = {}
        for A, B in zip(possible_remote_sites, prob_dist):
            prob_dict_remote_sites_given_local[str(A)] = B
        print("local sites are: " + str(local_sites))
        print("possible out of range sites are: " + str(possible_remote_sites))
        print("probability distribution is: " + str(prob_dist))
        for remote_sites in possible_remote_sites:
            if sites_invalid_for_agent_number(local_sites, remote_sites, num_agents):
                continue
            else:
                prob_remote_sites_given_local = prob_dict_remote_sites_given_local[str(remote_sites)]
                prob_sum = np.sum(list(prob_dict_remote_sites_given_local.values()))
                if(prob_sum != 1):
                    print("invalid distribution!")
                    break;
                prob_global_sites_given_motion = prob_remote_sites_given_local * prob_local_sites_given_motion_and_agent_num
                constructed_sites = np.concatenate((local_sites, remote_sites), axis=None)
                if str(constructed_sites) not in constructed_possible_sites:
                    constructed_possible_sites[str(constructed_sites)] = prob_global_sites_given_motion
                else: 
                    constructed_possible_sites[str(constructed_sites)] = constructed_possible_sites[str(constructed_sites)] + prob_global_sites_given_motion
                print("constructed global sites are: " + str(constructed_sites))
    return constructed_possible_sites
    
def get_mutual_information_given_motion_dict(global_params):
    # total_global_entropy = np.sum(np.array(list(weighted_global_entropy_given_motion.values())))
    total_macrostate_entropy = get_total_entropy_of_system_macrostate(global_params)
    energy_macrostate_given_motion_dict = get_global_entropy_of_macrostate_given_motion_dict(global_params)

    mutual_information_given_motion = dict()

    possible_motion = calculate_possible_sites(global_params["len_of_time_window"])
    for motion in possible_motion:
        macrostate_entropy_given_motion = energy_macrostate_given_motion_dict[str(motion)]
        mutual_info = total_macrostate_entropy - macrostate_entropy_given_motion
        mutual_information_given_motion[str(motion)] = mutual_info
    return mutual_information_given_motion


def get_observability_given_motion_dict(global_params):
    mutual_information_given_motion_dict = get_mutual_information_given_motion_dict(global_params)
    total_macrostate_entropy = get_total_entropy_of_system_macrostate(global_params)

    obsevability_given_motion_dict = dict()
    for motion in mutual_information_given_motion_dict.keys():
        mutual_information = mutual_information_given_motion_dict[motion]
        obsevability_given_motion_dict[motion] = mutual_information / total_macrostate_entropy
        # motion_array = convert_str_to_array(motion)
        # obsevability_given_motion_dict[motion] = calculate_observability(motion_array, total_macrostate_entropy, global_params)
    return obsevability_given_motion_dict
 

# def get_site_prob_dict_for_known_size_given_motion(motion, global_params):
#     possible_local_sites = calculate_possible_sites(global_params["influential_range"])
#     prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(motion, global_params)
#     # prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(no_recent_motion, global_params)
#     x_data = [str(local_sites) for local_sites in possible_local_sites]
#     y_data = [item.item() for item in prob_dist_local_sites_given_motion]
#     prob_dict_local_sites_given_motion = {}
#     for A, B in zip(x_data, y_data):
#         prob_dict_local_sites_given_motion[str(A)] = B

#     constructed_possible_sites = dict()
#     possible_num_of_agent = list(range(1,7))
#     for num_agents in possible_num_of_agent:
#         # prob_of_agent_num = 1 / len(possible_num_of_agent)
#         # prob_of_agent_num = calculate_prob_of_agent_num(num_agents, global_params['num_of_sites'])
#         prob_dict_local_sites_given_motion_and_agent_num = calculate_prob_dict_local_sites_given_motion_and_agent_num(prob_dict_local_sites_given_motion, num_agents, global_params)
#         # final_sum = np.sum(list(prob_dict_local_sites_given_motion_and_agent_num.values()))
#         for local_sites in possible_local_sites: 
#             # prob_of_local_agent_num = calculate_prob_of_local_agent_num(num_local_agents, possible_local_sites)
#             prob_local_sites_given_motion_and_agent_num = prob_dict_local_sites_given_motion_and_agent_num[str(local_sites)] #* prob_of_agent_num
#             # prob_local_sites_given_motion = prob_dict_local_sites_given_motion[str(local_sites)]
#             possible_remote_sites, prob_dist = calculate_prob_dist_remote_sites_given_local_sites_and_agent_num(local_sites, num_agents, global_params)
#             # make prob distribution uniform
#             # prob_dist[prob_dist >= 0] = 0.25
#             prob_dict_remote_sites_given_local = {}
#             for A, B in zip(possible_remote_sites, prob_dist):
#                 prob_dict_remote_sites_given_local[str(A)] = B

#             for remote_sites in possible_remote_sites:
#                 if sites_invalid_for_agent_number(local_sites, remote_sites, num_agents):
#                     continue
#                 else:
#                     prob_remote_sites_given_local = prob_dict_remote_sites_given_local[str(remote_sites)]
#                     prob_global_sites_given_motion = prob_remote_sites_given_local * prob_local_sites_given_motion_and_agent_num
#                     constructed_sites = np.concatenate((local_sites, remote_sites), axis=None)
#                     if str(constructed_sites) not in constructed_possible_sites:
#                         constructed_possible_sites[str(constructed_sites)] = prob_global_sites_given_motion
#                     else: 
#                         constructed_possible_sites[str(constructed_sites)] = constructed_possible_sites[str(constructed_sites)] + prob_global_sites_given_motion
#     return constructed_possible_sites
    