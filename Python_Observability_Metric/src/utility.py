# system libs
from distutils.command import clean
import math, itertools
from scipy import signal

# dependencies - numpy, matplotlib
import numpy as np 
from numpy.random import rand
from numpy import log, dot, e
np.set_printoptions(linewidth=np.inf)

def calculate_prob_dist(valued_array):
    """Calculate a probability distribution from the given
    array of values. (normalize the values so they all sum to 1)

    Args:
        valued_array (np.array): an array with numeric values which should be normalized

    Returns:
        np.array: 
            an array representing a probability distribution corresponding with the given values
    """
    prob_dist = valued_array / np.sum(valued_array[:])
    return prob_dist


def convert_array_to_unity_values(non_unity_array):
    # % localStates = (localStates==1).*1;
    func = lambda s: int(s==1)
    unity_array = np.array([func(xi) for xi in non_unity_array])
    return unity_array


def de2bi(d):
    """Converts a range of decimal values 'd' into a matrix of binary powers

    Args:
        d (np.array): an array consisting of a range of decimal values to be represented

    Returns:
        np.array: an array with rows consisting of binary representations of given decimal values
    """
    # calculate length of the rows in the returned matrix
    n = np.floor(np.sqrt(np.max(d))).astype(int) + 1;
    d = np.array(d)
    d = np.reshape(d, (1, -1))
    power = np.flipud(2**np.arange(n))

    g = np.zeros((np.shape(d)[1], n))

    for i, num in enumerate(d[0]):
        g[i] = num * np.ones((1,n))
    g = np.flip(g)
    b = np.floor((g%(2*power))/power)
    return b


def calculate_permutations_of_n_sites(num_of_sites):
    lst = list(map(list, itertools.product([0, 1], repeat=num_of_sites)))
    return np.array(lst)

def calculate_possible_sites(num_of_sites):
    """Calculates the permutations of all posible configurations of the sites 
        which fall within the influential range

    Args:
        num_of_sites (int): 
            an integer representing the range over which possible sites should be calculated

    Returns:
        [np.array]: 
            an num_of_sites x 2^(num_of_sites) size matrix representing every permution of possible 
            occucied (1) or unocupied (0) states for the range in consideration
    """
    # calculate the total number of possible state configurations that exist for the sites within the influential range
    # num_of_local_config = 2**num_of_sites;
    # create a temporary matrix representing every permutation of possible binary states
    # for a system with the given range
    # state_range = np.arange(0, (num_of_local_config ));
    # temp_state_mat = de2bi(state_range);
    temp_state_mat = calculate_permutations_of_n_sites(num_of_sites);
    # build a matrix representing every possible permutation of binary states for the length of the influential range
    # influential_range_vec = np.arange(1, influential_range + 1)
    # # influential_range_vec = np.arange(0, influential_range)
    # possible_local_sites = temp_state_mat[:, influential_range_vec];
    return temp_state_mat;


def get_prob_from_hamiltonion(H):
    """Given a value which represents a system Hamiltonion 
    calcuate a corresponding probability

    Args:
        H ([type]): a float rperesenting the hamiltonion (kinetic + potential energy)

    Returns:
        float: a float representing the corresponding probability
    """
    exp_fact = np.exp(-H)               # probability based on Hamiltonian
    P_tr = np.minimum(1, exp_fact);               # probability based on Hamiltonian
    return P_tr


# binom(n, k)
def nchoosek(n, k):
    """ Gives the number of possible ways that 'k' objects could be chosen
    out of a collection of 'n' total objects
    (this is irrespective of the order of the objects)

    Args:
        n (int): the total number of objects from which the choice is being made
        k (int): the number of objects being chosen

    Returns:
        int: the number of possible ways that 'k' objecs could be chosen out of 'n' total
    
    Examples:
    >>> nchoosek(10,6)
    210
    
    >>> nchoosek(1,1)
    1
    
    >>> nchoosek(5, 3)
    10
    """
    n = int(n)
    k = int(k)
    result = 0.0
    try:
        result = math.factorial(n) // math.factorial(k) // math.factorial(n - k)
    except:
        result = 0.0
    return result


## returns an array the length of the influential range
# which represents the interaction coefficients
def calculate_weighted_local_interactions(global_params):
    """Calculate what the strength of interactions should be for
    sites within the local range of an agent.

    Args:
        global_params (dict): the global simulation parameters

    Returns:
        np.array: an array representing the strength of interactions for
            agents within the influential range
    """
    interaction_coeff = global_params["interaction_coeff"]
    influential_range = global_params["influential_range"]
    influential_range_vec = np.arange(1, influential_range + 1)
    interaction_coeff_vec = np.zeros(influential_range + 1);
    # for each element in the influential range
    # populate the vector with a power of the hyperbolic tangent 'activation function'
    # this function takes input in any range and return a number in the range (-1, 1)
    for i in influential_range_vec:
        interaction_coeff_vec[i] = np.tanh(interaction_coeff)**i;
    interaction_coeff_vec = np.delete(interaction_coeff_vec, (0), axis=0)
    return interaction_coeff_vec

###
# calculates a gaussian distribution of the given length and steepness
# this can be used as a window to filter values when combined with some other function response
#
# takes:
# L - the length of the windo to calculate
# alpha - a parameter indicating how steep the curve should be
# returns:
# w - an array representing a gaussian distribution with a number of elements equal to L
###
def gausswin(L, alpha=2.5):
    """
    An N-point Gaussian window with alpha proportional to the
    reciprocal of the standard deviation.  The width of the window
    is inversely related to the value of alpha.  A larger value of
    alpha produces a more narrow window.
    Parameters
    ----------------------------
    L : int
    alpha : float
      Defaults to 2.5
    Returns
    ----------------------------
    Notes
    ----------------------------
    TODO: I am ignoring some corner cases, for example:
      #L - negative, error
      #L = 0
      #w => empty
      #L = 1
      #w = 1
    Equivalent of Matlab's gausswin function.
    """
    N = L - 1
    n = np.arange(0, N + 1) - N / 2
    w = np.exp(-(1 / 2) * (alpha * n / (N / 2)) ** 2)
    # return w
    return signal.windows.gaussian(L, alpha)


# This function helps determine the probability of the global states
def calculate_entropy_from_prob_dist(prob_dist):
    """Take a probability distribution and calculate an entropy from it
    the entropy here provides a measure of the average amount of information needed to represent 
    an event drawn from the probability distribution for a random variable

    Args:
        prob_dist (np.array): an array consisting of values in the (0-1) range
            representing a probability distribution

    Returns:
        float: a floating point number representing the corresponding entropy
            this is the average number of bits which is needed to represent an event
            from the distrubiton 
    """
    # ignore the 0 probability states
    non_zero_prob_dist = prob_dist[np.nonzero(prob_dist)]
    res = 0
    totalLength = non_zero_prob_dist.shape[0]
    for i in np.arange(0, totalLength):
        res = res - non_zero_prob_dist[i] * np.log2(non_zero_prob_dist[i])
    return res


def get_hamiltonion_for_site_from_local_sites(site, local_sites, global_params):
    external_field_coeff = global_params["external_field_coeff"]
    weighted_local_interactions = calculate_weighted_local_interactions(global_params).transpose()
    # calculate the internal field term by applying the weighted interaction coeffiecents to the local sites
    H_int = np.dot(weighted_local_interactions, local_sites)
    # External field term
    H_ext = external_field_coeff * site
    # compinee^{-\beta H_i} the internal and external hamiltonions to get the net result
    H = H_int - H_ext
    return H

def get_hamiltonion_for_system_sites(sites, global_params):
    H = 0
    for i, site in enumerate(sites):
        local_sites = get_local_sites_for_site(i, sites, global_params)
        site_ham = get_hamiltonion_for_site_from_local_sites(site, local_sites, global_params)
        H = H + site_ham
    return H

def get_local_sites_for_site(ind, system_sites, global_params):
    influential_range = global_params["influential_range"]
    num_of_sites = len(system_sites)
    # local_sites = np.zeros(influential_range)
    local_sites = np.ones(influential_range) * -1
    for i in np.arange(ind, ind + influential_range):
        if(i > num_of_sites - 1):
            local_sites[i-ind] = system_sites[i-num_of_sites]
        else:
            local_sites[i-ind] = system_sites[i]
    # local_sites = convert_array_to_unity_values(local_sites)
    return local_sites



def shift_binary_sites_to_negative(sites):
    """ shift sites consististing of 0,1 to -1,1  for use in hamiltonion calculattion

    Args:
        sites (np.array): the 0,1 sites to shift

    Returns:
        np.array: the -1,1 sites after shifting
    """
    shifted_sites = np.array(sites)*2 - 1
    return shifted_sites


def shift_negative_sites_to_binary(sites):
    """ shift sites consististing of -1,1 to 0,1  for use in analysis

    Args:
        sites (np.array): the -1,1 sites to shift

    Returns:
        np.array: the 0,1 sites after shifting
    """
    shifted_sites = np.array(sites)*0.5 + 1
    return shifted_sites

def clean_array_string(site_str_array):
    return str(site_str_array).strip(' []')


def convert_str_to_array(array_str):
    clean_str = clean_array_string(array_str)
    result = np.array([])
    for char in clean_str.split(' '):
        result = np.append(result, int(char))
    return result