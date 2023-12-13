# dependencies - numpy
from glob import glob
import numpy as np 

# local files
from probabilities import *
from utility import *

# %  Calculate the probability distributions by merging
# %  distributions within and beyond influential range
# %  P(\Sigma_complementary | \Sigma_m)* P(\Sigma_m| Y)
def calculate_probs_dist_global_sites_given_local_agents(num_local_agents, global_params):
    """[summary]
    TODO ???????????????????? need better understanding

    Args:
        num_of_local_agents (int): the number of agents which are found within the influential range
        global_params (dict): a dictionary containing the global simulation parameters

    Returns:
        np.array: an array containing probability values in the [0-1] range
        (this is a probability distrubution so its elements will sum to 1) 
    """
    num_of_agents = global_params["num_of_agents"]
    num_of_sites = global_params["num_of_sites"]

    # the extra '-1' here is due to the self exclusion of the agents from its own count
    num_of_agents_in_range = num_local_agents
    num_of_agents_out_of_range = num_of_agents - 1 - num_of_agents_in_range

    num_of_sites_in_range = global_params["influential_range"]
    num_of_sites_out_of_range = num_of_sites - 1 - num_of_sites_in_range

    prob_dist_out_of_range_agents = calculate_prob_dist_agents_out_of_range(num_of_agents_out_of_range, num_of_sites_out_of_range)
    prob_dist_in_range_agents = calculate_prob_dist_agents_in_sites(num_of_agents_in_range, num_of_sites_in_range)
    prob_global_dist_mat = np.tensordot(prob_dist_out_of_range_agents, prob_dist_in_range_agents, axes=0)
    prob_global_dist_vec = np.reshape(prob_global_dist_mat.transpose(),(1, -1),order='F')
    return prob_global_dist_vec


def estimate_global_entropy_given_num_of_local_agents(num_of_local_agents, global_params):
    """Calculate the overall global entropy of the system given a known number of agents present
    within the influential range
    this is a number which represents the total information needed to describe the possible system state
    given what is known about the local agents
    TODO ????????????????????? relies on function for which understanding is lacking currently

    Args:
        num_of_local_agents (int): the number of agents which are found within the influential range
        global_params (dict): a dictionary containing the global simulation parameters

    Returns:
        float: a floating point number representing the system entropy (ie: total information needed to describe the system)
    """
    probs_dist_global_sites = calculate_probs_dist_global_sites_given_local_agents(num_of_local_agents, global_params)
    entropy = calculate_entropy_from_prob_dist(probs_dist_global_sites) #H(\Sigma)
    return entropy


def estimate_global_entropy(global_params):
    """Calculate the overall global entropy of the system given the global parameters of the simulation
    this is a number which represents the total information needed to describe the possible system state

    Args:
        global_params (dict): a dictionary containing the global simulation parameters

    Returns:
        float: a floating point number representing the system entropy (ie: total information needed to describe the system)
    """
    # expand the global parameters
    num_of_agents = global_params["num_of_agents"]
    influential_range = global_params["influential_range"]
    # calculate the gloabl state entropy  H(\Simga)
    global_entropy = 0
    # for each possible number of agents occupying sites in the influential range 
    for num_of_local_agents in np.arange(0, min(num_of_agents, influential_range) + 1):
        # add an estimate of the entropy corresponding with this number of agents to the global entropy
        global_entropy = global_entropy + estimate_global_entropy_given_num_of_local_agents(num_of_local_agents, global_params)
    return global_entropy


def calculate_direct_global_entropy(global_params):
    """Calculate the overall global entropy of the system given the global parameters of the simulation
    this is a number which represents the total information needed to describe the possible system state

    Args:
        global_params (dict): a dictionary containing the global simulation parameters

    Returns:
        float: a floating point number representing the system entropy (ie: total information needed to describe the system)
    """
    # expand the global parameters
    num_of_agents = global_params["num_of_agents"]
    num_of_sites = global_params["num_of_sites"]
    
    prob_dist = calculate_prob_dist_agents_in_sites(num_of_agents, num_of_sites)
    entropy = calculate_entropy_from_prob_dist(prob_dist)

    return entropy

def calculate_probs_ff_given_possible_local_sites(possible_local_sites, global_params):
    """Calculate the probability for a free flow state given each possible configuration of sites within
    the influential range.

    Args:
        possible_local_sites (np.array): a multidimensional array with dimensions [influential_rang^2 x influential_range] 
            the rows of the matrix represent different configurations for the local sites (sites within influential range)
            and each column represents a site. the individual elements are valued (1,0) and represent whether the site
            is occupied or unoccupied
        global_params (dict): a dictionary containing the global simulation parameters

    Returns:
        np.array: an array of length [influential_rang^2] containing values in the range
            [0-1] which represent the probabilities of free-flow given each possible local site configuration
            (not a probability distribution, so wont sum to 1)
    """
    # calculate the interaction coefficient vector
    num_of_local_config = len(possible_local_sites)
    probs_ff_given_possible_local_sites = np.zeros((num_of_local_config, 1));
    # for i in np.arange(1, num_of_local_config):
    for i in np.arange(0, num_of_local_config):
        local_sites = possible_local_sites[i]
        probs_ff_given_possible_local_sites[i] = calculate_prob_ff_from_local_site_config(local_sites, global_params)
    return probs_ff_given_possible_local_sites

    
def calculate_prob_dict_ff_given_possible_local_sites(possible_local_sites, global_params):
    """Calculate the probability for a free flow state given each possible configuration of sites within
    the influential range.

    Args:
        possible_local_sites (np.array): a multidimensional array with dimensions [influential_rang^2 x influential_range] 
            the rows of the matrix represent different configurations for the local sites (sites within influential range)
            and each column represents a site. the individual elements are valued (1,0) and represent whether the site
            is occupied or unoccupied
        global_params (dict): a dictionary containing the global simulation parameters

    Returns:
        np.array: an array of length [influential_rang^2] containing values in the range
            [0-1] which represent the probabilities of free-flow given each possible local site configuration
            (not a probability distribution, so wont sum to 1)
    """
    # calculate the interaction coefficient vector
    num_of_local_config = len(possible_local_sites)
    # prob_dict_ff_given_local_sites = np.zeros((num_of_local_config, 1));
    prob_dict_ff_given_local_sites = dict();
    for local_sites in possible_local_sites:
        prob_dict_ff_given_local_sites[str(local_sites)] = calculate_prob_ff_from_local_site_config(local_sites, global_params)
    return prob_dict_ff_given_local_sites


def calculate_likelihood_local_sites_given_agent_motion(recent_agent_motion, global_params):
    """Calculate the likelihood for local sites to be in any configuration given the 
    mesured agent motion data

    Args:
        recent_agent_motion (np.array):
            an array representing a time window slice of an agents motion durring the simulation
            the array has elements of value (1,0) and a length of [len_of_time]
            the elements represent whether a given agent transitioned to the site in front of it
            after the given moment in time
        global_params (dict): a dictionary containing the global simulation parameters

    Returns:
        np.array: an array of length [influential_rang^2] containing floating point values which represent the estimated likelihood
            for any site configuration
    """
    influential_range = global_params["influential_range"]
    # calculate the probabilities of all possible local site configurations for the influential range
    possible_local_sites = calculate_possible_sites(influential_range)
    probs_of_possible_local_sites  = calculate_probs_of_possible_local_sites(influential_range)
    # calculate the probability of free flow from this
    prob_ff_from_recent_motion = calculate_prob_ff_from_recent_motion(recent_agent_motion.transpose())
    prob_jam_from_recent_motion = 1 - prob_ff_from_recent_motion
    # calculate the probabilitys for possible local sites given the motion data from the time window
    probs_ff_given_possible_local_sites = calculate_probs_ff_given_possible_local_sites(possible_local_sites, global_params);
    probs_jam_given_possible_local_sites = (1 - probs_ff_given_possible_local_sites)
    # calculate the probabilities for jam / free-flow given any of the possible local sites
    probs_ff_and_local_sites = probs_ff_given_possible_local_sites * probs_of_possible_local_sites
    probs_jam_and_local_sites = probs_jam_given_possible_local_sites * probs_of_possible_local_sites
    # from thes probabilities calculate the likelyhood for the local sites to be in any configuration
    probs_ff_and_local_sites_given_agent_motion = 0
    probs_jam_and_local_sites_given_agent_motion = 0
    # avoid dividing by 0 if any of the probabilities comes to that
    if(prob_ff_from_recent_motion != 0):
        probs_ff_and_local_sites_given_agent_motion = probs_ff_and_local_sites / prob_ff_from_recent_motion
    if(prob_jam_from_recent_motion != 0):
        probs_jam_and_local_sites_given_agent_motion = probs_jam_and_local_sites / prob_jam_from_recent_motion
    # the overall probability of the local site configuration is the combination of probabilities for free-flow and for jam
    likelihood_local_sites_given_agent_motion = probs_ff_and_local_sites_given_agent_motion + probs_jam_and_local_sites_given_agent_motion
    return likelihood_local_sites_given_agent_motion

# This function determines the likelihood of local states given measurement and conduct normalization
# derive_and_normalize_likelihood_from_bayesian
def calculate_prob_dist_local_sites_given_motion(recent_agent_motion, global_params):
    """ Calculate the probabilities for each possible local site configuration based 
    on what can be gleaned from the measured recent agent motion data

    Args:
        recent_agent_motion (np.array):
            an array representing a time window slice of an agents motion durring the simulation
            the array has elements of value (1,0) and a length of [len_of_time]
            the elements represent whether a given agent transitioned to the site in front of it
            after the given moment in time
        global_params (dict): a dictionary containing the global simulation parameters

    Returns:
        np.array: an array of length [influential_rang^2] containing a distribution of values in the range
            [0-1] which represent the probabilities of each possible local site configuration based on 
            the given agent motion data
    """
    likelihood_local_sites_given_agent_motion = calculate_likelihood_local_sites_given_agent_motion(recent_agent_motion, global_params)
    # normalise the likelihood arary to get a probability distribution for the result
    probs_local_sites_given_ff = calculate_prob_dist(likelihood_local_sites_given_agent_motion)
    return probs_local_sites_given_ff

def calculate_possible_out_of_range_sites_given_agent_num(possible_out_of_range_sites, global_params):
    num_agents = global_params["num_of_agents"]
    influential_range = global_params["influential_range"]
    possible_sites_given_agents = np.array([])
    for sites in list(possible_out_of_range_sites):
        if list(sites).count(1) >= (num_agents - influential_range):
            np.append(possible_sites_given_agents, sites)
    return possible_sites_given_agents
 
        
def calculate_possible_sites_given_agent_num(possible_local_sites, global_params):
    num_agents = global_params["num_of_agents"]
    influential_range = global_params["influential_range"]
    possible_sites_given_agents = np.array([])
    for sites in list(possible_local_sites):
        if list(sites).count(1) <= (num_agents - influential_range):
            np.append(possible_sites_given_agents, sites)
    return possible_sites_given_agents
    

def calculate_state_probs_from_recent_motion(recent_agent_motion):
    """Calculate the probabilites for both free-flow and jam states
    from the recent measured motion

    Args:
        recent_agent_motion (np.array):
            an array representing a time window slice of an agents motion durring the simulation
            the array has elements of value (1,0) and a length of [len_of_time]
            the elements represent whether a given agent transitioned to the site in front of it
            after the given moment in time

    Returns:
        float, float: two floating point numbers in the range (0-1) representing the 
            probabilities for both free-flow and jam states respectively given the 
            recent measured agent motion
    """
    prob_ff_from_recent_motion = calculate_prob_ff_from_recent_motion(recent_agent_motion.transpose())
    prob_jam_from_recent_motion = (1-prob_ff_from_recent_motion)
    return prob_ff_from_recent_motion, prob_jam_from_recent_motion

def calculate_prob_dist_remote_sites_given_local_sites_and_agent_num(local_sites, num_agents, global_params):
    # expand global parameters
    # num_agents = global_params["num_of_agents"]
    num_of_sites = global_params["num_of_sites"]
    influential_range = global_params["influential_range"]
    
    num_of_out_of_range_sites = num_of_sites - influential_range - 1
    num_local_agents = local_sites.tolist().count(1)
    if num_local_agents + num_of_out_of_range_sites + 1 < num_agents:
        num_of_out_of_range_agents = 0
    elif num_local_agents + 1 > num_agents:
        num_of_out_of_range_agents = 0
    else:
        num_of_out_of_range_agents = num_agents - num_local_agents - 1

    possible_out_of_range_sites = calculate_possible_sites(num_of_out_of_range_sites)
    
    # possible_sites_given_agents = calculate_possible_sites_given_agent_num(possible_out_of_range_sites, global_params)

    prob_dist = np.zeros(len(possible_out_of_range_sites))
    # for i in np.arange(0, num_of_out_of_range_sites):
    for i, out_of_range_sites in enumerate(possible_out_of_range_sites):
        if out_of_range_sites.tolist().count(1) == num_of_out_of_range_agents:
            num_of_ways_to_fit_agents_in_sites = nchoosek(num_of_out_of_range_sites, num_of_out_of_range_agents)
            prob_dist[i] = 1 / num_of_ways_to_fit_agents_in_sites
    if np.sum(prob_dist) != 1.0:
        print("TEST")
    return possible_out_of_range_sites, prob_dist

def calculate_probs_out_of_range_agents_given_states(local_sites_data, global_params):
    """Calculate probabilities for the out of range agents ....
    TODO need to figure this out

    Args:
        local_sites_data (np.array): 
            a multidimensional array with a length [2^influential_range] ie: the number of possible 
            local site configurations. each element contains both an array of length [influential array]
            which contains the corresponding specific site configuration and probability calculated for
            this configuration from the motion data
        global_params (dict): a dictionary containing the global simulation parameters

    Returns:
        np.array, np.array: two arrays representing a probability distribution for
            agents out of the influential range correpsponding with both free-flow
            and jam states
    """
    # expand global parameters
    num_of_agents = global_params["num_of_agents"]
    num_of_sites = global_params["num_of_sites"]
    influential_range = global_params["influential_range"]
    # expand the local site data
    local_sites = local_sites_data[0]
    prob_local_sites = local_sites_data[1]
    # get the number of agents occupying these local sites
    num_local_agents = local_sites.tolist().count(1)
    # calculate the numbor of agents and sites which are outside of range
    num_of_agents_out_of_range = num_of_agents - 1 - num_local_agents
    num_of_sites_out_of_range = num_of_sites - 1 - influential_range
    # calculate the probability distribution which corresponds with this
    prob_dist_out_of_range_agents_given_local_sites = calculate_prob_dist_agents_out_of_range(num_of_agents_out_of_range, num_of_sites_out_of_range)
    probs_out_of_range_agents_given_ff = prob_dist_out_of_range_agents_given_local_sites * prob_local_sites;
    probs_out_of_reach_agents_given_jam = 1-probs_out_of_range_agents_given_ff
    return probs_out_of_range_agents_given_ff, probs_out_of_reach_agents_given_jam
    

def calculate_entropy_for_local_sites_given_agent_motion(recent_agent_motion, global_params, local_sites_data):
    """Calculate the global system entropy given a time slice of the motion data for an agent and a specific
    site configuration.

    Args:
        recent_agent_motion (np.array):
            an array representing a time window slice of an agents motion durring the simulation
            the array has elements of value (1,0) and a length of [len_of_time]
            the elements represent whether a given agent transitioned to the site in front of it
            after the given moment in time
        global_params (dict): the global simulation parameters
        local_sites_data (np.array): 
            a multidimensional array with a length [2^influential_range] ie: the number of possible 
            local site configurations. each element contains both an array of length [influential array]
            which contains the corresponding specific site configuration and probability calculated for
            this configuration from the motion data

    Returns:
        float: a floating point number representing the amount of information needed to describe the possible
            system based on what can be known from the motion history of the agent in question
    """
    # calculate the probabilities for both the free-flow and jam states from the measurement motion data
    prob_ff_from_recent_motion, prob_jam_from_recent_motion  = calculate_state_probs_from_recent_motion(recent_agent_motion)
    # calculate the probability distrubitons for the out of range eagents to be in these states
    probs_out_of_range_agents_given_ff, probs_out_of_range_agents_given_jam = calculate_probs_out_of_range_agents_given_states(local_sites_data, global_params)
    # calculate the entropies for both the jam and free-flow states from these probability distrubutions
    entropy_for_out_of_range_agents_given_ff = calculate_entropy_from_prob_dist(probs_out_of_range_agents_given_ff)
    entropy_for_out_of_range_agents_given_jam = calculate_entropy_from_prob_dist(probs_out_of_range_agents_given_jam)
    # TODO figure this out
    # MATLAB - % 
    # ProbabilityOfFreeFlowFromData*obj.entropyCalculation(Probability_SigmaGlobal_Given_Y) + 
    # (1-ProbabilityOfFreeFlowFromData)*obj.entropyCalculation(1-Probability_SigmaGlobal_Given_Y)
    # scale these entropy values by the probabilities of their dependent states given the motion measurement
    local_sites_entropy = prob_ff_from_recent_motion * entropy_for_out_of_range_agents_given_ff
    other_sites_entropy = prob_jam_from_recent_motion * entropy_for_out_of_range_agents_given_jam
    # calculate the total entropy by combining these
    total_entropy = local_sites_entropy  + other_sites_entropy
    # return the result
    return total_entropy

def calculate_system_entropy_given_measurement(recent_agent_motion, global_params):
    """Calculate the sum of the conditional entropy for the global system
    given a certain time slice of an agents motion during the simulation
    the conditional entropy represents the information needed to describe the outcome 
    of one unknown random variable (the sites of the system) given the known values 
    of another random variable (the measured ff values from the time slice)

    Args:
        recent_agent_motion (np.array):
            an array representing a time window slice of an agents motion durring the simulation
            the array has elements of value (1,0) and a length of [len_of_time]
            the elements represent whether a given agent transitioned to the site in front of it
            after the given moment in time
        global_params (dict): a dictionary containing the global simulation parameters

    Returns:
        float: a floating point number representing the ammount of information 
            needed to describe the possible system given the measured motion data
    """
    influential_range = global_params["influential_range"]
    possible_local_sites = calculate_possible_sites(influential_range)
    prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(recent_agent_motion, global_params);
    possible_local_sites_data = np.array(list(zip(possible_local_sites,prob_dist_local_sites_given_motion)))

    sum_of_entropy_for_sites_given_ff = 0
    # iterate all possible configuration
    # for index in np.arange(0, num_of_local_sites_config):
    for local_sites_data in possible_local_sites_data:
        # Calculation of conditional entropy   H(\Sigma|Y)
        entropy = calculate_entropy_for_local_sites_given_agent_motion(recent_agent_motion, global_params, local_sites_data)
        sum_of_entropy_for_sites_given_ff  = sum_of_entropy_for_sites_given_ff + entropy
    return sum_of_entropy_for_sites_given_ff


def calculate_mutual_information(recent_agent_motion, global_entropy, global_params):
    """ Calculate the mutual information value which correpsonds
    with the given agent motion measurement data captured for a time window

    Args:
        recent_agent_motion (np.array):
            an array representing a time window slice of an agents motion durring the simulation
            the array has elements of value (1,0) and a length of [len_of_time]
            the elements represent whether a given agent transitioned to the site in front of it
            after the given moment in time
        global_params (dict): the global simulation parameters

    Returns:
        float: a number representing how observable the global state of the system is given the measured agent motion data
    """
    # calculate entropy of the system given the measured agent motion 
    # (this represents the information needed to describe the system given the measurements)
    sum_of_entropy_for_sites_given_ff = calculate_system_entropy_given_measurement(recent_agent_motion, global_params)
    # definition of mutual information
    # the mutual information represents the information which we know from the measurd agent motion 
    # this is the difference between the overall information needed for the system and the information needed given our measurement
    mutual_information = global_entropy - sum_of_entropy_for_sites_given_ff 
    # probabilityMeasurementDistribution = np.array([prob_ff_from_recent_motion, 1-prob_ff_from_recent_motion])
    # entropyOfMeasurement = calculate_entropy_from_probs(probabilityMeasurementDistribution)   #H(Y)
    # measurement_efficiency_metric = mutual_information/entropyOfMeasurement
    return mutual_information

# This function helps determine the probability of the global states
def calculate_observability(recent_agent_motion, global_entropy, global_params):
    """ Calculate the observability metric value which correpsonds
    with the given agent motion measurement data captured for a time window

    Args:
        recent_agent_motion (np.array):
            an array representing a time window slice of an agents motion durring the simulation
            the array has elements of value (1,0) and a length of [len_of_time]
            the elements represent whether a given agent transitioned to the site in front of it
            after the given moment in time
        global_params (dict): the global simulation parameters

    Returns:
        float: a number representing how observable the global state of the system is given the measured agent motion data
    """
    # calculate entropy of the system given the measured agent motion 
    # (this represents the information needed to describe the system given the measurements)
    sum_of_entropy_for_sites_given_ff = calculate_system_entropy_given_measurement(recent_agent_motion, global_params)
    # definition of mutual information
    # the mutual information represents the information which we know from the measurd agent motion 
    # this is the difference between the overall information needed for the system and the information needed given our measurement
    mutual_information = global_entropy - sum_of_entropy_for_sites_given_ff 
    # probabilityMeasurementDistribution = np.array([prob_ff_from_recent_motion, 1-prob_ff_from_recent_motion])
    # entropyOfMeasurement = calculate_entropy_from_probs(probabilityMeasurementDistribution)   #H(Y)
    # measurement_efficiency_metric = mutual_information/entropyOfMeasurement

    # the observability is then the ration of mutual_information to overall information
    # it has value 1 if mutual info = overall info, and value 0 if mutual info = 0
    observability = mutual_information/global_entropy

    return observability


def get_observability_metric_series_for_agent(agent_motion_through_time, global_params, index_of_agent):
    """ Create a time series containing observability measurements for each point in time for
    which these measurements can be attained with the time window specified in the simulation parameters

    Args:
        agent_motion_through_time (np.array): 
            an array representing agent motion over the course of the simulation
            the array has elements of value (1,0) and dimensions [num_of_agents x len_of_time]
            the elements represent whether a given agent transitioned to the site in front of it
            after the given moment in time
        global_params (dict): the global simulation parameters
        index_of_agent (int): the index of the agent under consideration

    Returns:
        np.array: an array with length [len_of_time - len_of_time_window] whose elements 
            represent the observability measures for corresponding points in time
    """
    len_of_time = global_params["len_of_time"]
    len_of_time_window = global_params["len_of_time_window"]
    
    seriesLenth = len_of_time - len_of_time_window
    observability_series = np.zeros((seriesLenth, 1))
    t0 = len_of_time_window
    # global_entropy = estimate_global_entropy(global_params)
    global_entropy = calculate_direct_global_entropy(global_params)
    for t in np.arange(t0, t0+seriesLenth):
    # t0 = 0
    # calculate the overall entropy of the system 
    # (this represents the information needed to descripbe the systems state)
    # for t in np.arange(t0, seriesLenth):
        # acquire this agents motion data for the given time window
        # recent_agent_motion = agent_motion_through_time[index_of_agent, t:t+len_of_time_window]
        recent_agent_motion = agent_motion_through_time[index_of_agent, t-len_of_time_window:t]
        # finally calculate the observability
        observability = calculate_observability(recent_agent_motion, global_entropy, global_params)
        # add this to the series
        time_ind = t-t0
        observability_series[time_ind] = observability
    return observability_series