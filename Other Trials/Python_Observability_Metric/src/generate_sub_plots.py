# system libs
import sys
import math

from random import Random
random = Random(1337)

# dependencies - numpy, matplotlib
import numpy as np 
from numpy.random import rand
from numpy import NaN, log, dot, e
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage,AnnotationBbox

from plotly.subplots import make_subplots
import plotly.graph_objects as go

# local files
from utility import *
from simulation import *
from analysis import *
from probabilities import *
from plot_utilities import *
from example_analysis import *

np.set_printoptions(threshold=sys.maxsize)

# Traffic Model prameters
# and Parameters for information quantification
global_params = {
    "num_of_sites": 6,
    "num_of_agents": 3,
    "interaction_coeff": 5.0,
    "external_field_coeff": 0.8,
    "len_of_time": 64,
    "influential_range": 3,
    "len_of_time_window": 4
}


len_of_time_window = global_params['len_of_time_window']
full_motion_data = np.ones(len_of_time_window).astype(float)
empty_motion_data = np.zeros(len_of_time_window).astype(float)

# folder = "figures_time" + str(global_params["len_of_time"]) + "_sites" + str(global_params["num_of_sites"]) + "_m" + str(global_params["influential_range"]) +  "_J" + str(global_params["interaction_coeff"]) + "_F" + str(global_params["external_field_coeff"])

#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################


fig, axs = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle(r"$p(\mathcal{Y}^k=y_f) = \sum_{i=k-m}^k \sigma_0^i g(i),  \quad$ where $\:  \: g(i) = \frac{\sqrt{2}}{\sqrt{\pi s^2}}exp(-\frac{i^2}{2s^2})$", fontsize=16)

possible_motion = calculate_possible_sites(len_of_time_window)
x_data = []
y_data = []
for motion_history in possible_motion: 
    res = calculate_prob_ff_from_recent_motion(motion_history)
    y_data.append(res)
    x_data.append(str(motion_history))
series_data = zip(x_data, y_data)
series_list = list(series_data)
series_list.sort(key=lambda x:x[1])
unzipped_data = list(zip(*series_list))
axs.bar(unzipped_data[0], 1.0, color="firebrick", label="probability of jam")
rects = axs.bar(unzipped_data[0], unzipped_data[1], label="probability of ff")
axs.tick_params('x', labelrotation=90)
axs.set_title('probability ff given observed motion')
axs.set_xlabel('Observed Motion: Y \n (oldest at bottom)')
axs.set_ylabel('Probability of Free-flow')
axs.legend()
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axs.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]), 2)), ha='center', va='bottom')
plt.savefig(f"./Python_Observability_Metric/figures/prob_ff_given_motion.pdf")


#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################


fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle(r"$p(\Sigma_m = {\sigma}_{m_i})$", fontsize=16)

# calculate the probabilities of all possible local site configurations for the influential range
possible_local_sites = calculate_possible_sites(global_params['influential_range'])
probs_of_possible_local_sites  = calculate_probs_of_possible_local_sites(global_params['influential_range'])
x_data = [str(local_sites) for local_sites in possible_local_sites]
y_data = [item[0] for item in probs_of_possible_local_sites]
series_data = zip(x_data, y_data)
series_list = list(series_data)
# series_list.sort(key=lambda x:x[1])
series_list.sort(key=lambda x: int(str(x[0]).strip(' []').replace(',','').replace(' ',''),2))
unzipped_data = list(zip(*series_list))
rects = axes.bar(unzipped_data[0], unzipped_data[1])
axes.tick_params('x', labelrotation=90)
axes.set_title('Probability of Local Site Configuration')
axes.set_xlabel('Local site configuration \n ')
axes.set_ylabel('Probability of local sites')
axes.legend()
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]), 2)), ha='center', va='bottom')

# for i, site_str in enumerate(x_data):
sorted_x = list(zip(*series_list))[0]
for i, site_str in enumerate(sorted_x):
    site_config_label = site_str[::-1].strip(' []').replace(',','').replace(' ','')
    offset_image(i, site_config_label, axes)
plt.savefig("./Python_Observability_Metric/figures/prob_local_config.pdf")


#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################


fig, axes = plt.subplots(1, 1, constrained_layout=True)
# fig.suptitle(r"$p(\mathbf{Y}) = \sum_i p(\mathbf{Y}|\Sigma_m = \underline{\sigma}_{m_i})p(\Sigma_m = \underline{\sigma}_{m_i})$", fontsize=16)
# fig.suptitle(r"$p(\mathbf{Y}) = \sum_i p(\mathbf{Y}|\Sigma_m = {\sigma}_{m_i})p(\Sigma_m = {\sigma}_{m_i})$", fontsize=16)
fig.suptitle(r"$p(\mathcal{Y}^k=y_f | \Sigma_m = {\sigma}_{m_i}))$", fontsize=16)
    
possible_local_sites = calculate_possible_sites(global_params["influential_range"])
x_data = []
y_data = []
for local_sites in possible_local_sites: 
    res = calculate_prob_ff_from_local_site_config(local_sites, global_params)
    y_data.append(res)
    x_data.append(str(local_sites))
series_data = zip(x_data, y_data)
series_list = list(series_data)
series_list.sort(key=lambda x:x[0], reverse=True)
series_list.sort(key=lambda x:x[1])
unzipped_data = list(zip(*series_list))
axes.bar(unzipped_data[0], 1.0, color='firebrick', label="probability of jam")
rects = axes.bar(unzipped_data[0], unzipped_data[1], label="probability of ff")
axes.tick_params('x', labelrotation=90)
axes.set_title('probability ff given local sites')
axes.set_xlabel('Local Site configuration \n ')
axes.set_ylabel('Probability of Free-flow')
axes.legend()
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),2)), ha='center', va='bottom')

sorted_x = list(zip(*series_list))[0]
for i, site_str in enumerate(sorted_x):
    site_config_label = site_str[::-1].strip(' []').replace(',','').replace(' ','')
    offset_image(i, site_config_label, axes)
plt.savefig("./Python_Observability_Metric/figures/prob_ff_given_local_config.pdf")


#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################


fig, axes = plt.subplots(2, 1, constrained_layout=True)
fig.suptitle(r"$\mathcal{L}(\Sigma_m = \sigma_{m_i}|\mathbf{Y}) = \frac{p(\mathbf{Y}|\Sigma_m = \sigma_{m_i})p(\Sigma_m = \sigma_{m_i})}{p(\mathbf{Y})}$", fontsize=16)

possible_local_sites = calculate_possible_sites(global_params["influential_range"])
full_recent_motion = np.ones(global_params["len_of_time_window"]).astype(float)
res = calculate_likelihood_local_sites_given_agent_motion(full_recent_motion, global_params)
x_data = [str(local_sites) for local_sites in possible_local_sites]
y_data = [item.item() for item in res]
# for local_sites in possible_local_sites: 
#     y_data.append(res)
#     x_data.append(str(local_sites))
series_data = zip(x_data, y_data)
series_list = list(series_data)
series_list.sort(key=lambda x:x[0], reverse=True)
series_list.sort(key=lambda x:x[1])
unzipped_data = list(zip(*series_list))
rects = axes[0].bar(unzipped_data[0], unzipped_data[1])
axes[0].tick_params('x', labelrotation=90)
axes[0].set_title(r"likeliehood of local sites given full recent motion")
axes[0].set_xlabel('Local Site configuration \n ')
axes[0].set_ylabel('Likeliehood')
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axes[0].text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')

# for i, site_str in enumerate(x_data):
sorted_x = list(zip(*series_list))[0]
for i, site_str in enumerate(sorted_x):
    site_config_label = site_str[::-1].strip(' []').replace(',','').replace(' ','')
    offset_image(i, site_config_label, axes[0])

possible_local_sites = calculate_possible_sites(global_params["influential_range"])
full_recent_motion = np.zeros(global_params["len_of_time_window"]).astype(float)
res = calculate_likelihood_local_sites_given_agent_motion(full_recent_motion, global_params)
x_data = [str(local_sites) for local_sites in possible_local_sites]
y_data = [item.item() for item in res]
# for local_sites in possible_local_sites: 
#     y_data.append(res)
#     x_data.append(str(local_sites))
series_data = zip(x_data, y_data)
series_list = list(series_data)
series_list.sort(key=lambda x:x[0], reverse=True)
series_list.sort(key=lambda x:x[1], reverse=True)
unzipped_data = list(zip(*series_list))
rects = axes[1].bar(unzipped_data[0], unzipped_data[1])
axes[1].tick_params('x', labelrotation=90)
axes[1].set_title('likeliehood of local sites given no recent motion')
axes[1].set_xlabel('Local Site configuration \n ')
axes[1].set_ylabel('Likeliehood')
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axes[1].text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')

fig.tight_layout(rect=[0, 0.03, 1, 0.95])

# for i, site_str in enumerate(x_data):
sorted_x = list(zip(*series_list))[0]
for i, site_str in enumerate(sorted_x):
    site_config_label = site_str[::-1].strip(' []').replace(',','').replace(' ','')
    offset_image(i, site_config_label, axes[1])

plt.savefig("./Python_Observability_Metric/figures/likelihood_local_sites_given_motion.pdf")


#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################


fig, axes = plt.subplots(2, 1, constrained_layout=True)
# fig.suptitle(r"$p(\Sigma_m = \underline{\sigma}_i|\mathbf{Y}) = \frac{\mathcal{L}(\Sigma_m = \underline{\sigma}_i|\mathbf{Y})}{\sum_i  \mathcal{L}(\Sigma_m = \underline{\sigma}_i|\mathbf{Y})}$", fontsize=16)
fig.suptitle(r"$p(\Sigma_m = {\sigma}_i|\mathbf{Y}) = \frac{\mathcal{L}(\Sigma_m = {\sigma}_i|\mathbf{Y})}{\sum_i  \mathcal{L}(\Sigma_m = {\sigma}_i|\mathbf{Y})}$", fontsize=16)

empty_recent_motion = np.ones(global_params["len_of_time_window"]).astype(float)
res = calculate_prob_dist_local_sites_given_motion(empty_recent_motion, global_params)
x_data = [str(local_sites) for local_sites in possible_local_sites]
y_data = [item.item() for item in res]
# for local_sites in possible_local_sites: 
#     y_data.append(res)
#     x_data.append(str(local_sites))
series_data = zip(x_data, y_data)
series_list = list(series_data)
series_list.sort(key=lambda x:x[0], reverse=True)
series_list.sort(key=lambda x:x[1])
unzipped_data = list(zip(*series_list))
rects = axes[0].bar(unzipped_data[0], unzipped_data[1])
axes[0].tick_params('x', labelrotation=90)
axes[0].set_title('probability of local sites given full recent motion')
axes[0].set_xlabel('Local Site configuration \n ')
axes[0].set_ylabel('probability')
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axes[0].text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')

# for i, site_str in enumerate(x_data):
sorted_x = list(zip(*series_list))[0]
for i, site_str in enumerate(sorted_x):
    site_config_label = site_str[::-1].strip(' []').replace(',','').replace(' ','')
    offset_image(i, site_config_label, axes[1])
    site_config_label = site_str[::-1].strip(' []').replace(',','').replace(' ','')
    offset_image(i, site_config_label, axes[0])



empty_recent_motion = np.zeros(global_params["len_of_time_window"]).astype(float)
res = calculate_prob_dist_local_sites_given_motion(empty_recent_motion, global_params)
x_data = [str(local_sites) for local_sites in possible_local_sites]
y_data = [item.item() for item in res]
# for local_sites in possible_local_sites: 
#     y_data.append(res)
#     x_data.append(str(local_sites))
series_data = zip(x_data, y_data)
series_list = list(series_data)
series_list.sort(key=lambda x:x[0], reverse=True)
series_list.sort(key=lambda x:x[1], reverse=True)
unzipped_data = list(zip(*series_list))
rects = axes[1].bar(unzipped_data[0], unzipped_data[1])
axes[1].tick_params('x', labelrotation=90)
axes[1].set_title('probability of local sites given no recent motion')
axes[1].set_xlabel('Local Site configuration \n ')
axes[1].set_ylabel('probability')
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axes[1].text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')

fig.tight_layout()

# for i, site_str in enumerate(x_data):
sorted_x = list(zip(*series_list))[0]
for i, site_str in enumerate(sorted_x):
    site_config_label = site_str[::-1].strip(' []').replace(',','').replace(' ','')
    offset_image(i, site_config_label, axes[1])
    site_config_label = site_str[::-1].strip(' []').replace(',','').replace(' ','')
    offset_image(i, site_config_label, axes[1])

plt.savefig("./Python_Observability_Metric/figures/prob_local_sites_given_motion.pdf")


#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################


# fig, axes = plt.subplots(1,1, constrained_layout=True)
# fig.suptitle(r"$H(\Sigma |\mathbf{Y}) = \sum_{y\in \mathbf{Y}} p(y) H(\Sigma|\mathbf{Y} = y)$", fontsize=16)


# possible_motion = calculate_possible_sites(len_of_time_window)
# x_data = []
# y_data = []
# for motion_history in possible_motion: 
#     res = calculate_system_entropy_given_measurement(motion_history, global_params)
#     y_data.append(res)
#     x_data.append(str(motion_history))
# series_data = zip(x_data, y_data)
# series_list = list(series_data)
# series_list.sort(key=lambda x:x[1])
# unzipped_data = list(zip(*series_list))
# rects = axes.bar(unzipped_data[0], unzipped_data[1])
# axes.tick_params('x', labelrotation=90)
# axes.set_title('System entropy given observed motion')
# axes.set_xlabel('Observed Motion: Y  \n (oldest at bottom)')
# axes.set_ylabel('System Entropy')
# # add text to the graph
# for i, rect in enumerate(rects):
#     height = rect.get_height()
#     axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),2)), ha='center', va='bottom')
# plt.savefig("./Python_Observability_Metric/figures/entropy.pdf")


# plt.savefig("./Python_Observability_Metric/figures/original_entropy_given_motion.pdf")


#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################


# fig, axes = plt.subplots(1, 1, constrained_layout=True)
# fig.suptitle(r"$I(\Sigma;\mathbf{Y}) = H(\Sigma)-H(\Sigma|\mathbf{Y})$", fontsize=16)

# possible_motion = calculate_possible_sites(len_of_time_window)
# x_data = []
# y_data = []
# global_entropy = calculate_direct_global_entropy(global_params)
# # global_entropy = estimate_global_entropy(global_params)
# for motion_history in possible_motion: 
#     res = calculate_mutual_information(motion_history, global_entropy, global_params)
#     y_data.append(res)
#     x_data.append(str(motion_history))
# series_data = zip(x_data, y_data)
# series_list = list(series_data)
# series_list.sort(key=lambda x:x[1])
# unzipped_data = list(zip(*series_list))
# rects = axes.bar(unzipped_data[0], unzipped_data[1])
# axes.tick_params('x', labelrotation=90)
# axes.set_title('Mutual Information of system given observed motion')
# axes.set_xlabel('Observed Motion: Y \n (oldest at bottom)')
# axes.set_ylabel('Mutual Information')
# # add text to the graph
# for i, rect in enumerate(rects):
#     height = rect.get_height()
#     axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),2)), ha='center', va='bottom')

# plt.savefig("./Python_Observability_Metric/figures/original_mutual_information_given_motion.pdf")


#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################


# sup_title = r"$\rho(\Sigma|\mathbf{Y}) = \frac{I(\Sigma;\mathbf{Y})}{H(\Sigma)} $"
# graph = MatplotGraph(1, 1, sup_title)
# possible_motion = calculate_possible_sites(len_of_time_window)
# global_entropy = calculate_direct_global_entropy(global_params)

# def generate_y(motion_history):
#     res = calculate_observability(motion_history, global_entropy, global_params)
#     return motion_history, res
# title = 'Observability of system given observed motion'
# x_label = 'Observed Motion: Y \n (oldest at bottom)'
# y_label = 'Observability'
# graph.plot_subgraph(possible_motion, generate_y, x_label, y_label, title, 0, False, 'alphanumeric', 0)
# graph.save_to_file("./Python_Observability_Metric/figures/original_observability_given_motion.pdf")

#######################################################################
#######################################################################
#######################################################################


# sup_title = r"$H(\Sigma_m | \mathbf{Y}) = H(\Sigma_m | y_f)P(y_f | \mathbf{Y}) + H(\Sigma_m | y_j)P(y_j | \mathbf{Y})$"
# graph = MatplotGraph(2, 1, sup_title)
# possible_local_sites = calculate_possible_sites(global_params["influential_range"])

# prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(full_motion_data, global_params);
# possible_local_sites_data = np.array(list(zip(possible_local_sites,prob_dist_local_sites_given_motion)))
# def generate_y_for_full(local_sites_data):
#     local_sites = local_sites_data[0]
#     res = calculate_entropy_for_local_sites_given_agent_motion(full_motion_data, global_params, local_sites_data)
#     return local_sites, res
# title = 'entropy of local sites given full recent motion'
# x_label = 'Local Site configuration \n '
# y_label = 'Entropy'
# graph.plot_subgraph(possible_local_sites_data, generate_y_for_full, x_label, y_label, title, 0, True, 'binary', 0)

# prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(empty_motion_data, global_params);
# possible_local_sites_data = np.array(list(zip(possible_local_sites,prob_dist_local_sites_given_motion)))
# def generate_y_for_empty(local_sites_data):
#     local_sites = local_sites_data[0]
#     res = calculate_entropy_for_local_sites_given_agent_motion(empty_motion_data, global_params, local_sites_data)
#     return local_sites, res
# title = 'entropy of local sites given no recent motion'
# x_label = 'Local Site configuration \n '
# y_label = 'Entropy'
# graph.plot_subgraph(possible_local_sites_data, generate_y_for_empty, x_label, y_label, title, 1, True, 'binary', 0)
# graph.save_to_file("./Python_Observability_Metric/figures/entropy_local_sites_given_motion.pdf")


####################################################################################################################
####################################################################################################################
####################################################################################################################


sup_title = r"$\mathcal{H} | \Sigma = \sum_i (-K\sum_j \sigma_i \sigma_j - B \sigma_i )$"
graph = MatplotGraph(1, 1, sup_title)
possible_sites = calculate_possible_sites(global_params["num_of_sites"])

def generate_y(possible_sites):
    # shifted_possible_sites = shift_binary_sites_to_negative(possible_sites)
    res = get_hamiltonion_for_system_sites(possible_sites, global_params)
    cleaned_sites = clean_array_string(possible_sites)
    return cleaned_sites, res
title = 'Effective Hamiltonion for system sites'
x_label = 'Global site configuration'
y_label = 'Hamiltonion'
graph.plot_subgraph(possible_sites, generate_y, x_label, y_label, title, 0, True, 'alphanumeric', 1)
graph.set_figure_size(6, 20)
graph.save_to_file("./Python_Observability_Metric/figures/hamiltonion_given_all_sites.pdf")
 

####################################################################################################################
####################################################################################################################
####################################################################################################################


# fig, axes = plt.subplots(1, 1, constrained_layout=True)
# # fig.suptitle(r"$p(\Sigma_m = \underline{\sigma}_i|\mathbf{Y}) = \frac{\mathcal{L}(\Sigma_m = \underline{\sigma}_i|\mathbf{Y})}{\sum_i  \mathcal{L}(\Sigma_m = \underline{\sigma}_i|\mathbf{Y})}$", fontsize=16)
# # fig.suptitle(r"$p(\Sigma | \mathbf{Y}) = p(\Sigma_{m}, \Sigma_{m}^*|\mathbf{Y}) = p(\Sigma_{m}^* | \Sigma_{m}) \cdot p(\Sigma_{m}| \mathbf{Y})$", fontsize=16)
# fig.suptitle(r"$p(\Sigma | \mathbf{Y}) = p(\Sigma_{m}, \Sigma_{m}^*|\mathbf{Y}) = p(\Sigma_{m}^* | \Sigma_{m}) \cdot p(\Sigma_{m}| \mathbf{Y})$ given 3 agents and full motion", fontsize=16)


# full_recent_motion = np.ones(global_params["len_of_time_window"]).astype(float)

# num_agents = 3
# constructed_possible_sites = get_site_prob_dict_given_motion_and_agent_num(full_recent_motion, num_agents, global_params)

# x_data = constructed_possible_sites.keys()
# y_data = constructed_possible_sites.values()

# plt.rcParams.update({'font.size': 6})
# series_data = zip(x_data, y_data)
# series_list = list(series_data)
# series_list.sort(key=lambda x:x[0], reverse=True)
# series_list.sort(key=lambda x:x[1])

# unzipped_data = list(zip(*series_list))
# rects = axes.bar(unzipped_data[0], unzipped_data[1])
# axes.tick_params('x', labelrotation configuration')xend 3otions
# axes.set_title('probability of site configuration given full agent motion and 3 agents')
# axes.set_xlabel('site configuration \n ')
# axes.set_ylabel('probability')
# # add text to the graph
# for i, rect in enumerate(rects):
#     height = rect.get_height()
    
# # for i, site_str in enumerate(x_data):
# sorted_x = list(zip(*series_list))[0]
# for i, site_str in enumerate(sorted_x):
#     site_config_label = site_str[::-1].strip(' []').replace(',','').replace(' ','')
#     offset_image(i, site_config_label, axes)
    
# # spacing = 0.400
# # fig.subplots_adjust(bottom=spacing)

# t(y_data))th(12)ig, axes y_data))
# print("sum of global = pl possibilities given full motion and 3 agents is: " + str(prob_sum))

# fig.set_figheight(8)
# fig.set_figwidth(12)
# # fig.set_sizetinches(5, 4)
# # plt.tight_layout()
# plt.savefig("./Python_Observability_Metric/figures/.sub_pll_sites_given_full_motion_3_agents.pdf")

#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################

fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle(r"$p(\Sigma | \mathbf{Y}) = p(\Sigma_{m}, \Sigma_{m}^*|\mathbf{Y}) = p(\Sigma_{m}^* | \Sigma_{m}) \cdot p(\Sigma_{m}| \mathbf{Y})$", fontsize=16)

possible_local_sites = calculate_possible_sites(global_params["influential_range"])
full_recent_motion = np.ones(global_params["len_of_time_window"]).astype(float)
no_recent_motion = np.zeros(global_params["len_of_time_window"]).astype(float)

constructed_possible_sites = get_site_prob_given_motion_dict(full_recent_motion, global_params)

x_data = constructed_possible_sites.keys()
site_probabilities = constructed_possible_sites.values()

plt.rcParams.update({'font.size': 6})
series_data = zip(x_data, site_probabilities)
series_list = list(series_data)
series_list.sort(key=lambda x:x[0], reverse=True)
series_list.sort(key=lambda x:x[1])

# series_list.sort(key=lambda x: int(str(x[0]).strip(' []').replace(',','').replace(' ',''),2))
unzipped_data = list(zip(*series_list))
rects = axes.bar(unzipped_data[0], unzipped_data[1])
axes.tick_params('x', labelrotation=90)
axes.set_title('probability of all sites given full agent motion')
# axes.set_title('probability of global sites given no agent motion')
axes.set_xlabel('site configuration \n ')
axes.set_ylabel('probability')
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    # axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')
prob_sum = np.sum(list(site_probabilities))
p_sum = np.sum(list(site_probabilities))
print("sum of global site possibilities given full agent motion comes to: " + str(prob_sum))
# 
# for i, site_str in enumerate(x_data):
sorted_x = list(zip(*series_list))[0]
for i, site_str in enumerate(sorted_x):
    # site_config_label = site_str[::-1].strip(' []').replace(',','').replace(' ','')
    site_config_label = site_str.strip(' []').replace(',','').replace(' ','')
    offset_image(i, site_config_label, axes)


fig.set_figheight(8)
fig.set_figwidth(12)
plt.savefig("./Python_Observability_Metric/figures/prob_global_sites_given_full_motion.pdf")


#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################


fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle(r"$p(\Sigma | \mathbf{Y}) = p(\Sigma_{m}, \Sigma_{m}^*|\mathbf{Y}) = p(\Sigma_{m}^* | \Sigma_{m}) \cdot p(\Sigma_{m}| \mathbf{Y})$", fontsize=16)

half_recent_motion = np.zeros(global_params["len_of_time_window"]).astype(float)
half_recent_motion[::2] += 1
# prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(full_recent_motion, global_params)

prob_system_microstate_given_motion_dict = get_site_prob_given_motion_dict(half_recent_motion, global_params)

x_data = prob_system_microstate_given_motion_dict.keys()
site_probabilities = prob_system_microstate_given_motion_dict.values()

# plt.rcParams.update({'font.size': 6})
series_data = zip(x_data, site_probabilities)
series_list = list(series_data)
# series_list.sort(key=lambda x:x[0], reverse=True)

series_list.sort(key=lambda x: int(str(x[0]).strip(' []').replace(',','').replace(' ',''),2), reverse=True)
series_list.sort(key=lambda x:x[1], reverse=False)
unzipped_data = list(zip(*series_list))
rects = axes.bar(unzipped_data[0], unzipped_data[1])
axes.tick_params('x', labelrotation=90)
axes.set_title('probability of global sites given half agent motion')
axes.set_xlabel('site configuration \n ')
axes.set_ylabel('probability')
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    # axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')
prob_sum = np.sum(list(site_probabilities))
print("sum of global site possibilities given half agent motion comes to: " + str(prob_sum))
# 
# for i, site_str in enumerate(x_data):
sorted_x = list(zip(*series_list))[0]
for i, site_str in enumerate(sorted_x):
    # site_config_label = site_str[::-1].strip(' []').replace(',','').replace(' ','')
    site_config_label = site_str[::].strip(' []').replace(',','').replace(' ','')
    offset_image(i, site_config_label, axes)


fig.set_figheight(8)
fig.set_figwidth(12)
plt.savefig("./Python_Observability_Metric/figures/prob_global_sites_given_half_motion.pdf")


#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################


fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle(r"$p(\Sigma | \mathbf{Y}) = p(\Sigma_{m}, \Sigma_{m}^*|\mathbf{Y}) = p(\Sigma_{m}^* | \Sigma_{m}) \cdot p(\Sigma_{m}| \mathbf{Y})$", fontsize=16)

full_recent_motion = np.ones(global_params["len_of_time_window"]).astype(float)
no_recent_motion = np.zeros(global_params["len_of_time_window"]).astype(float)
# prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(full_recent_motion, global_params)

prob_system_microstate_given_motion_dict = get_site_prob_given_motion_dict(no_recent_motion, global_params)

x_data = prob_system_microstate_given_motion_dict.keys()
site_probabilities = prob_system_microstate_given_motion_dict.values()
# 
# plt.rcParams.update({'font.size': 6})
series_data = zip(x_data, site_probabilities)
series_list = list(series_data)
series_list.sort(key=lambda x:x[0], reverse=True)
series_list.sort(key=lambda x:x[1], reverse=True)

# series_list.sort(key=lambda x: int(str(x[0]).strip(' []').replace(',','').replace(' ',''),2))
unzipped_data = list(zip(*series_list))
rects = axes.bar(unzipped_data[0], unzipped_data[1])
axes.tick_params('x', labelrotation=90)
axes.set_title('probability of global sites given no agent motion')
axes.set_xlabel('site configuration \n ')
axes.set_ylabel('probability')
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    # axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')
prob_sum = np.sum(list(site_probabilities))
print("sum of global site possibilities given no agent motion comes to: " + str(prob_sum))
# 
# for i, site_str in enumerate(x_data):
sorted_x = list(zip(*series_list))[0]
for i, site_str in enumerate(sorted_x):
    # site_config_label = site_str[::-1].strip(' []').replace(',','').replace(' ','')
    site_config_label = site_str[::].strip(' []').replace(',','').replace(' ','')
    offset_image(i, site_config_label, axes)


fig.set_figheight(8)
fig.set_figwidth(12)
plt.savefig("./Python_Observability_Metric/figures/prob_global_sites_given_no_motion.pdf")



#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################

fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle(r"$p(E_{\Sigma} | \mathbf{Y}) = p(\Sigma_{m}, \Sigma_{m}^*|\mathbf{Y}) = p(\Sigma_{m}^* | \Sigma_{m}) \cdot p(\Sigma_{m}| \mathbf{Y})$", fontsize=16)

possible_local_sites = calculate_possible_sites(global_params["influential_range"])
no_recent_motion = np.zeros(global_params["len_of_time_window"]).astype(float)

possible_energy_probs = get_energy_prob_dict_given_motion(no_recent_motion, global_params)
prob_dist = np.array(list(possible_energy_probs.values()))

x_data = possible_energy_probs.keys()
site_probabilities = possible_energy_probs.values()

plt.rcParams.update({'font.size': 6})
series_data = zip(x_data, site_probabilities)
series_list = list(series_data)
# series_list.sort(key=lambda x:x[0])
series_list.sort(key=lambda x: float(str(x[0])))
unzipped_data = list(zip(*series_list))
rects = axes.bar(unzipped_data[0], unzipped_data[1])
axes.tick_params('x', labelrotation=90)
axes.set_title('probability of system energy given no agent motion')
# axes.set_title('probability of global sites given no agent motion')
axes.set_xlabel('energy level \n ')
axes.set_ylabel('probability')
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')
prob_sum = np.sum(list(site_probabilities))

fig.set_figheight(8)
fig.set_figwidth(12)
plt.savefig("./Python_Observability_Metric/figures/prob_system_energy_given_no_motion.pdf")



#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################

fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle(r"$p(E_{\Sigma} | \mathbf{Y}) = p(\Sigma_{m}, \Sigma_{m}^*|\mathbf{Y}) = p(\Sigma_{m}^* | \Sigma_{m}) \cdot p(\Sigma_{m}| \mathbf{Y})$", fontsize=16)

possible_local_sites = calculate_possible_sites(global_params["influential_range"])
full_recent_motion = np.ones(global_params["len_of_time_window"]).astype(float)
no_recent_motion = np.zeros(global_params["len_of_time_window"]).astype(float)
# 
# constructed_possible_sites = get_site_prob_given_motion_dict(full_recent_motion, global_params)
possible_energy_probs = get_energy_prob_dict_given_motion(full_recent_motion, global_params)
prob_dist = np.array(list(possible_energy_probs.values()))

x_data = possible_energy_probs.keys()
site_probabilities = possible_energy_probs.values()

plt.rcParams.update({'font.size': 6})
series_data = zip(x_data, site_probabilities)
series_list = list(series_data)
# series_list.sort(key=lambda x:x[0], reverse=True)
series_list.sort(key=lambda x: float(str(x[0])))
unzipped_data = list(zip(*series_list))
rects = axes.bar(unzipped_data[0], unzipped_data[1])
axes.tick_params('x', labelrotation=90)
axes.set_title('probability of system energy given full agent motion')
# axes.set_title('probability of global sites given no agent motion')
axes.set_xlabel('energy level \n ')
axes.set_ylabel('probability')
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')
prob_sum = np.sum(list(site_probabilities))

fig.set_figheight(8)
fig.set_figwidth(12)
plt.savefig("./Python_Observability_Metric/figures/prob_system_energy_given_full_motion.pdf")


 
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################
#------------------------------------------------------------------------------------------------------------------ entropy of system microstate

fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle(r"$H(\Sigma | \mathbf{Y}) = \sum_{i} p(\Sigma_i | \mathbf{Y}) \log_2 p(\Sigma_i | \mathbf{Y}) $", fontsize=16)

global_entropy_given_motion = get_global_entropy_of_microstate_given_motion_dict(global_params)

plt.rcParams.update({'font.size': 6})
series_data = zip(list(global_entropy_given_motion.keys()), list(global_entropy_given_motion.values()))
series_list = list(series_data)
# series_list.sort(key=lambda x:x[0])
unzipped_data = list(zip(*series_list))
rects = axes.bar(unzipped_data[0], unzipped_data[1])
axes.tick_params('x', labelrotation=90)
axes.set_title("System entropy given observed motion")
# axes.set_title('probability of global sites given no agent motion')
axes.set_xlabel(r"observed motion $\mathbf{Y}$")
axes.set_ylabel(r"entopy of system $H(\Sigma | \mathbf{Y})$")
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')

plt.savefig("./Python_Observability_Metric/figures/entropy_given_motion.pdf")

 
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################
#------------------------------------------------------------------------------------------------------------------ entropy of system energy

fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle(r"$H(E_{\Sigma} | \mathbf{Y}) = \sum_{i} p(E_{\Sigma_i} | \mathbf{Y}) \log_2 p(E_{\Sigma_i} | \mathbf{Y}) $", fontsize=16)

global_entropy_macrostate_given_motion_dict = get_global_entropy_of_macrostate_given_motion_dict(global_params)

plt.rcParams.update({'font.size': 6})
series_data = zip(list(global_entropy_macrostate_given_motion_dict.keys()), list(global_entropy_macrostate_given_motion_dict.values()))
series_list = list(series_data)
# series_list.sort(key=lambda x:x[0])
unzipped_data = list(zip(*series_list))
rects = axes.bar(unzipped_data[0], unzipped_data[1])
axes.tick_params('x', labelrotation=90)
axes.set_title("System macrostate entropy given observed motion")
# axes.set_title('probability of global sites given no agent motion')
axes.set_xlabel(r"observed motion $\mathbf{Y}$")
axes.set_ylabel(r"entropy of system macrostate $H(E_{\Sigma} | \mathbf{Y})$")
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')

plt.savefig("./Python_Observability_Metric/figures/entropy_energy_given_motion.pdf")


 
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################
#------------------------------------------------------------------------------------------------------------------ entropy of system energy

fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle(r"$\Delta H(E_{\Sigma} | \mathbf{Y}) = H(E_{\Sigma} | \mathbf{Y}) - H(E_{\Sigma} | \mathbf{Y_0}) $", fontsize=16)

global_entropy_macrostate_given_motion_dict = get_global_entropy_of_macrostate_diff_given_motion_dict(global_params)

plt.rcParams.update({'font.size': 6})
series_data = zip(list(global_entropy_macrostate_given_motion_dict.keys()), list(global_entropy_macrostate_given_motion_dict.values()))
series_list = list(series_data)
series_list.sort(key=lambda x:x[0])
# series_list.sort(key=lambda x: int(str(x[0]).strip(' []').replace(',','').replace(' ',''),2))
unzipped_data = list(zip(*series_list))
rects = axes.bar(unzipped_data[0], unzipped_data[1])
axes.tick_params('x', labelrotation=90)
axes.set_title("Relative system macrostate entropy given observed motion")
# axes.set_title('probability of global sites given no agent motion')
axes.set_xlabel(r"observed motion $\mathbf{Y}$")
axes.set_ylabel(r"entropy gain of system macrostate $H(E_{\Sigma} | \mathbf{Y})$")
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')

plt.savefig("./Python_Observability_Metric/figures/entropy_energy_diff_given_motion.pdf")
 
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################
# ----------------------------------------------------------------------------------------- mutual information

fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle(r"$I(E_{\Sigma};\mathbf{Y}) = H(E_{\Sigma})-H(E_{\Sigma}|\mathbf{Y})$", fontsize=16)

mutual_information_given_motion = get_mutual_information_given_motion_dict(global_params)

plt.rcParams.update({'font.size': 6})
series_data = zip(list(mutual_information_given_motion.keys()), list(mutual_information_given_motion.values()))
series_list = list(series_data)
# series_list.sort(key=lambda x:x[0])
series_list.sort(key=lambda x: int(str(x[0]).strip(' []').replace(',','').replace(' ',''),2))
unzipped_data = list(zip(*series_list))
rects = axes.bar(unzipped_data[0], unzipped_data[1])
axes.tick_params('x', labelrotation=90)
axes.set_title("Mutual information given observed motion")
# axes.set_title('probability of global sites given no agent motion')
axes.set_xlabel(r"observed motion $\mathbf{Y}$")
axes.set_ylabel(r"mutual information of system $I(\Sigma;\mathbf{Y}) $")
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')

plt.savefig("./Python_Observability_Metric/figures/new_mutual_information_given_motion.pdf")

 
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################
# -----------------------------------------------------------------------------------------  observability

fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle(r"$\rho(E_{\Sigma} | \mathbf{Y}) = \frac{I(\Sigma;\mathbf{Y})}{H(E_{\Sigma})}$", fontsize=16)

observability_given_motion_dict = get_observability_given_motion_dict(global_params)

plt.rcParams.update({'font.size': 6})
series_data = zip(list(observability_given_motion_dict.keys()), list(observability_given_motion_dict.values()))
series_list = list(series_data)
# series_list.sort(key=lambda x:x[0])
series_list.sort(key=lambda x: int(str(x[0]).strip(' []').replace(',','').replace(' ',''),2))
unzipped_data = list(zip(*series_list))
rects = axes.bar(unzipped_data[0], unzipped_data[1])
axes.tick_params('x', labelrotation=90)
axes.set_title("Observability of system state given observed motion")
# axes.set_title('probability of global sites given no agent motion')
axes.set_xlabel(r"observed motion $\mathbf{Y}$")
axes.set_ylabel(r"odservability of system $\rho(E_{\Sigma} | \mathbf{Y}) $")
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')

plt.savefig("./Python_Observability_Metric/figures/new_observability_given_motion.pdf")

 
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################

fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle(r"$p(n) = \frac{\mathit{number \ of \ sites \ with \ n \ agents}}{\mathit{total \ number \ of \ sites}}$", fontsize=16)
x_data = []
y_data = []
# global_entropy = estimate_global_entropy(global_params)
possible_num_of_agent = list(range(1,7))
num_sites = global_params['num_of_sites']
for num_agents in possible_num_of_agent:
    prob_of_agent_num = calculate_prob_of_agent_num(num_agents, num_sites)
    y_data.append(prob_of_agent_num)
    x_data.append(num_agents)
series_data = zip(x_data, y_data)
series_list = list(series_data)
series_list.sort(key=lambda x:x[1])
unzipped_data = list(zip(*series_list))
rects = axes.bar(unzipped_data[0], unzipped_data[1])
axes.tick_params('x', labelrotation=90)
axes.set_title('Probability unknown system configuration contains N Agents')
axes.set_xlabel('number of agents')
axes.set_ylabel('probability of system with agents')
# add text to the graph
for i, rect in enumerate(rects):
    height = rect.get_height()
    axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),2)), ha='center', va='bottom')

plt.savefig("./Python_Observability_Metric/figures/prob_sytem_with_n_agents.pdf")


#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################
# -----------------------------------------------------------------------------------------  observability for different site size and observation window

fig, axes = plt.subplots(3, 2, constrained_layout=True)
fig.suptitle("Observability of system at different scales", fontsize=16)
fig.supxlabel("number of system sites")
fig.supylabel("size of observation window")

# Traffic Model prameters
# and Parameters for information quantification
global_params = {
    "num_of_sites": 4,
    "num_of_agents": 3,
    "interaction_coeff": 5.0,
    "external_field_coeff": 0.8,
    "len_of_time": 64,
    "influential_range": 3,
    "len_of_time_window": 2
}
observability_given_motion_dict = get_observability_given_motion_dict(global_params)
plt.rcParams.update({'font.size': 6})
series_data = zip(list(observability_given_motion_dict.keys()), list(observability_given_motion_dict.values()))
series_list = list(series_data)
# series_list.sort(key=lambda x:x[0])
series_list.sort(key=lambda x: int(str(x[0]).strip(' []').replace(',','').replace(' ',''),2))
unzipped_data = list(zip(*series_list))
rects = axes[0,0].bar(unzipped_data[0], unzipped_data[1])
axes[0,0].tick_params('x', labelrotation=90)
axes[0,0].set_title("Observability of system with 4 sites")
# axes[0,0].set_title('probability of global sites given no agent motion')
axes[0,0].set_xlabel(r"observed motion with time window 2")
axes[0,0].set_ylim([0,1])
axes[0,0].set_ylim([0,1])



global_params = {
    "num_of_sites": 12,
    "num_of_agents": 3,
    "interaction_coeff": 5.0,
    "external_field_coeff": 0.8,
    "len_of_time": 64,
    "influential_range": 3,
    "len_of_time_window": 2
}
observability_given_motion_dict = get_observability_given_motion_dict(global_params)
plt.rcParams.update({'font.size': 6})
series_data = zip(list(observability_given_motion_dict.keys()), list(observability_given_motion_dict.values()))
series_list = list(series_data)
# series_list.sort(key=lambda x:x[0])
series_list.sort(key=lambda x: int(str(x[0]).strip(' []').replace(',','').replace(' ',''),2))
unzipped_data = list(zip(*series_list))
rects = axes[0,1].bar(unzipped_data[0], unzipped_data[1])
axes[0,1].tick_params('x', labelrotation=90)
axes[0,1].set_title("Observability of system with 12 sites")
axes[0,1].set_xlabel(r"observed motion with time window 2")
axes[0,1].set_ylim([0,1])
axes[0,1].set_ylim([0,1])



global_params = {
    "num_of_sites": 4,
    "num_of_agents": 3,
    "interaction_coeff": 5.0,
    "external_field_coeff": 0.8,
    "len_of_time": 64,
    "influential_range": 3,
    "len_of_time_window": 4
}
observability_given_motion_dict = get_observability_given_motion_dict(global_params)
plt.rcParams.update({'font.size': 6})
series_data = zip(list(observability_given_motion_dict.keys()), list(observability_given_motion_dict.values()))
series_list = list(series_data)
# series_list.sort(key=lambda x:x[0])
series_list.sort(key=lambda x: int(str(x[0]).strip(' []').replace(',','').replace(' ',''),2))
unzipped_data = list(zip(*series_list))
rects = axes[1,0].bar(unzipped_data[0], unzipped_data[1])
axes[1,0].tick_params('x', labelrotation=90)
axes[1,0].set_title("Observability of system with 4 site5")
axes[1,0].set_xlabel(r"observed motion with time window 4")
axes[1,0].set_ylim([0,1])
axes[1,0].set_ylim([0,1])



global_params = {
    "num_of_sites": 12,
    "num_of_agents": 3,
    "interaction_coeff": 5.0,
    "external_field_coeff": 0.8,
    "len_of_time": 64,
    "influential_range": 3,
    "len_of_time_window": 4
}
observability_given_motion_dict = get_observability_given_motion_dict(global_params)
plt.rcParams.update({'font.size': 6})
series_data = zip(list(observability_given_motion_dict.keys()), list(observability_given_motion_dict.values()))
series_list = list(series_data)
# series_list.sort(key=lambda x:x[0])
series_list.sort(key=lambda x: int(str(x[0]).strip(' []').replace(',','').replace(' ',''),2))
unzipped_data = list(zip(*series_list))
rects = axes[1,1].bar(unzipped_data[0], unzipped_data[1])
axes[1,1].tick_params('x', labelrotation=90)
axes[1,1].set_title("Observability of system with 12 sites")
axes[1,1].set_xlabel(r"observed motion with time window 4")
axes[1,1].set_ylim([0,1])
axes[1,1].set_ylim([0,1])


global_params = {
    "num_of_sites": 4,
    "num_of_agents": 3,
    "interaction_coeff": 5.0,
    "external_field_coeff": 0.8,
    "len_of_time": 64,
    "influential_range": 3,
    "len_of_time_window": 5
}
observability_given_motion_dict = get_observability_given_motion_dict(global_params)
plt.rcParams.update({'font.size': 6})
series_data = zip(list(observability_given_motion_dict.keys()), list(observability_given_motion_dict.values()))
series_list = list(series_data)
# series_list.sort(key=lambda x:x[0])
series_list.sort(key=lambda x: int(str(x[0]).strip(' []').replace(',','').replace(' ',''),2))
unzipped_data = list(zip(*series_list))
rects = axes[2,0].bar(unzipped_data[0], unzipped_data[1])
axes[2,0].tick_params('x', labelrotation=90)
axes[2,0].set_title("Observability of system with 4 sites")
axes[2,0].set_xlabel(r"observed motion with time window 5")
axes[2,0].set_ylabel(r"odservability")
axes[2,0].set_ylim([0,1])



global_params = {
    "num_of_sites": 12,
    "num_of_agents": 3,
    "interaction_coeff": 5.0,
    "external_field_coeff": 0.8,
    "len_of_time": 64,
    "influential_range": 3,
    "len_of_time_window": 5
}
observability_given_motion_dict = get_observability_given_motion_dict(global_params)
plt.rcParams.update({'font.size': 6})
series_data = zip(list(observability_given_motion_dict.keys()), list(observability_given_motion_dict.values()))
series_list = list(series_data)
# series_list.sort(key=lambda x:x[0])
series_list.sort(key=lambda x: int(str(x[0]).strip(' []').replace(',','').replace(' ',''),2))
unzipped_data = list(zip(*series_list))
rects = axes[2,1].bar(unzipped_data[0], unzipped_data[1])
axes[2,1].tick_params('x', labelrotation=90)
axes[2,1].set_title("Observability of system with 12 sites")
axes[2,1].set_xlabel(r"observed motion with time window 5")
axes[2,1].set_ylabel(r"odservability")
axes[2,1].set_ylim([0, 1])



plt.savefig("./Python_Observability_Metric/figures/observability_scaling.pdf")