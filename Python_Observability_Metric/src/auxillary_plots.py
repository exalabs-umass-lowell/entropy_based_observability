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

np.set_printoptions(threshold=sys.maxsize)

# Traffic Model prameters
# and Parameters for information quantification
global_params = {
    "num_of_sites": 6,
    "num_of_agents": 3,
    "interaction_coeff": 4.0,
    "external_field_coeff": 1.5,
    "len_of_time": 64,
    "influential_range": 3,
    "len_of_time_window": 4
}


len_of_time_window = 4
full_motion_data = np.ones(len_of_time_window).astype(float)
empty_motion_data = np.zeros(len_of_time_window).astype(float)



#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################
# -------------------------------------------------------------------------------------------------------------------------------------------

# fig, axes = plt.subplots(1, 1, constrained_layout=True)
# # fig.suptitle(r"$\rho(E_{\Sigma}|\mathbf{Y}) = \frac{I(E_{\Sigma};\mathbf{Y})}{H(E_{\Sigma})} $", fontsize=16)
# x_data = []
# y_data = []
# # global_entropy = estimate_global_entropy(global_params)
# unique_macrostates =  np.unique(list(site_probabilities))
# for macrostate in unique_macrostates: 
#     prob_of_macrostate = list(site_probabilities).count(macrostate) / len(site_probabilities)
#     y_data.append(prob_of_macrostate)
#     x_data.append(str(round(macrostate, 4)))
# series_data = zip(x_data, y_data)
# series_list = list(series_data)
# series_list.sort(key=lambda x:x[1])
# unzipped_data = list(zip(*series_list))
# rects = axes.bar(unzipped_data[0], unzipped_data[1])
# axes.tick_params('x', labelrotation=90)
# axes.set_title('Probability of Macrostate')
# axes.set_xlabel('Macrostate')
# axes.set_ylabel('Probability of Macrostate')
# # add text to the graph
# for i, rect in enumerate(rects):
#     height = rect.get_height()
#     axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),2)), ha='center', va='bottom')

# print("sum of macrostate probabilities comes to: " + str(np.sum(list(y_data))))

# plt.savefig("prob_macrostate.svg")


#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################


tab_fig, tab_axes = plt.subplots(1, 1, constrained_layout=True)
tab_fig.suptitle(r"$p(n | \Sigma_{m}) = \frac{\mathit{number \ of \ site \ configurations \ with \ n \ agents}}{\mathit{total \ number \ of \ site \ configurations}}$", fontsize=16)
rows = []
cell_text = []
# global_entropy = estimate_global_entropy(global_params)
possible_local_sites = calculate_possible_sites(global_params['influential_range'])
possible_num_of_agent = list(range(1,7))
num_sites = global_params['num_of_sites'] 
for local_sites in possible_local_sites:
    row_cols = []
    for num_agents in possible_num_of_agent:
        prob_of_agent_num_given_local_sites = calculate_prob_agent_num_given_local_sites(num_agents, num_sites, local_sites)
        row_cols.append(prob_of_agent_num_given_local_sites)
    rows.append(str(local_sites))
    cell_text.append(row_cols)
columns = possible_num_of_agent
# Add a table at the bottom of the tab_axes
tab_axes.axis('tight')
tab_axes.axis('off')
tab_axes.set_title('Probability of agent number given local sites')
tab_axes.set_xlabel('number of agents')
tab_axes.set_ylabel('local site configuration')
the_table = tab_axes.table(cellText=cell_text,rowLabels=rows, colLabels=columns, loc='center')


# possible_local_sites = calculate_possible_sites(global_params["influential_range"])
# prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(motion, global_params)
# x_data = [str(local_sites) for local_sites in possible_local_sites]
# y_data = [item.item() for item in prob_dist_local_sites_given_motion]
# prob_dict_local_sites_given_motion = {}
# for A, B in zip(x_data, y_data):
#     prob_dict_local_sites_given_motion[str(A)] = B

# constructed_possible_sites = dict()
# possible_num_of_agent = list(range(0,global_params["num_of_agents"] + 1))
# for num_agents in possible_num_of_agent:
#     # prob_of_agent_num = 1 / len(possible_num_of_agent)
#     # prob_of_agent_num = calculate_prob_of_agent_num(num_agents, global_params['num_of_sites'])
#     # prob_dict_local_sites_given_motion_and_agent_num = calculate_prob_dict_local_sites_given_motion_and_agent_num(prob_dict_local_sites_given_motion, num_agents, global_params)
#     # final_sum = np.sum(list(prob_dict_local_sites_given_motion_and_agent_num.values()))
#     for local_sites in possible_local_sites: 
#         num_local_agents = list(local_sites).count(1)
#         # prob_of_local_agent_num = calculate_prob_of_local_agent_num(num_local_agents, possible_local_sites)
#         # prob_local_sites_given_motion_and_agent_num = prob_dict_local_sites_given_motion_and_agent_num[str(local_sites)] #* prob_of_agent_num
#         prob_local_sites_given_motion = prob_dict_local_sites_given_motion[str(local_sites)]
#         possible_remote_sites, prob_dist = calculate_prob_dist_remote_sites_given_local_sites_and_agent_num(local_sites, num_agents, global_params)
#         # prob_dist[prob_dist >= 0] = 0.25
#         prob_dict_remote_sites_given_local = {}
#         for A, B in zip(possible_remote_sites, prob_dist):
#             prob_dict_remote_sites_given_local[str(A)] = B

#         for remote_sites in possible_remote_sites:
#             num_remote_sites = len(possible_remote_sites)
#             num_remote_agents = list(remote_sites).count(1)
#             if num_local_agents + num_remote_agents + 1 != num_agents:
#                 continue
#             prob_remote_sites_given_local = prob_dict_remote_sites_given_local[str(remote_sites)]
#             prob_global_sites_given_motion = prob_remote_sites_given_local * prob_local_sites_given_motion
#             constructed_sites = np.concatenate((local_sites, remote_sites), axis=None)
#             if str(constructed_sites) not in constructed_possible_sites:
#                 constructed_possible_sites[str(constructed_sites)] = prob_global_sites_given_motion
#             else: 
#                 constructed_possible_sites[str(constructed_sites)] = constructed_possible_sites[str(constructed_sites)] + prob_global_sites_given_motion
# return constructed_possible_sites


#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################

# series_data = zip(rows, columns)
# series_list = list(series_data)
# series_list.sort(key=lambda x:x[1])
# unzipped_data = list(zip(*series_list))
# rects = axes.bar(unzipped_data[0], unzipped_data[1])
# add text to the graph




# fig14, fig14_axes = plt.subplots(1, 1, constrained_layout=True)
# fig14.suptitle(r"$p(\Sigma | \mathbf{Y}) = p(\Sigma_{m}, \Sigma_{m}^*|\mathbf{Y}) = p(\Sigma_{m}^* | \Sigma_{m}) \cdot p(\Sigma_{m}| \mathbf{Y})$", fontsize=16)

# possible_num_of_agent = list(range(0,global_params['influential_range']))
# for num_local_agents in possible_num_of_agent:
#     prob_dist = calculate_probs_dist_global_sites_given_local_agents(num_local_agents, global_params)

# x_data = possible_num_of_agent
# y_data = prob_dist

# series_data = zip(x_data, y_data)
# series_list = list(series_data)
# series_list.sort(key=lambda x:x[1])
# unzipped_data = list(zip(*series_list))
# rects = fig14_axes.bar(unzipped_data[0], unzipped_data[1])
# fig14_axes.tick_params('x', labelrotation=90)
# fig14_axes.set_title('probability of global sites given full agent motion')
# fig14_axes.set_xlabel('site configuration \n ')
# fig14_axes.set_ylabel('probability')
# # add text to the graph
# for i, rect in enumerate(rects):
#     height = rect.get_height()
#     fig14_axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),4)), ha='center', va='bottom')

# prob_sum = np.sum(list(site_probabilities))
# print("sum of global site possibilities given full agent motion comes to: " + str(prob_sum))






#######################################################################
#######################################################################
#######################################################################

sup_title = r"$\mathcal{H}_i = -K\sum_j \sigma_i \sigma_j - B \sigma_i $"
graph = MatplotGraph(2, 1, sup_title)
possible_local_sites = calculate_possible_sites(global_params["influential_range"])

def generate_y_for_full(local_sites):
    site = 1
    # shifted_sites = shift_binary_sites_to_negative(local_sites)
    res = get_hamiltonion_for_site_from_local_sites(site, local_sites, global_params)
    return local_sites, res
title = 'Effective Hamiltonian for occupied site ($\sigma_i = 1$)'
x_label = 'Local site configuration $\sigma_j$\n'
y_label = 'Hamiltonion'
graph.plot_subgraph(possible_local_sites, generate_y_for_full, x_label, y_label, title, 0, True, 'alphanumeric', 1)

def generate_y_for_empty(local_sites):
    site = 0
    # shifted_sites = shift_binary_sites_to_negative(local_sites)
    res = get_hamiltonion_for_site_from_local_sites(site, local_sites, global_params)
    return local_sites, res
title = 'Effective Hamiltonian for empty site ($\sigma_i = -1$)'
x_label = 'Local site configuration $\sigma_j$\n'
y_label = 'Hamiltonion'
graph.plot_subgraph(possible_local_sites, generate_y_for_empty, x_label, y_label, title, 1, True, 'alphanumeric', 1)

graph.save_to_file("./figures/hamiltonion_given_occupied_or_empty.pdf")






# # fig5, fig5_axes = plt.subplots(2, 1, constrained_layout=True)

# # possible_local_sites = calculate_possible_sites(global_params["influential_range"])
# # prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(full_motion_data, global_params);
# # possible_local_sites_data = np.array(list(zip(possible_local_sites,prob_dist_local_sites_given_motion)))
# # x_data = []
# # y_data = []
# # for local_sites_data in possible_local_sites_data:
# #     local_sites = local_sites_data[0]
# #     x_data.append(str(local_sites))
# #     num_of_agents = np.count_nonzero(local_sites)
# #     y_data.append(calculate_agents_energy_level_count(num_of_agents, global_params['num_of_sites']))
# # series_data = zip(x_data, y_data)
# # series_list = list(series_data)
# # # series_list.sort(key=lambda x:x[1])
# # unzipped_data = list(zip(*series_list))
# # rects = fig5_axes[0].bar(unzipped_data[0], unzipped_data[1])
# # fig5_axes[0].tick_params('x', labelrotation=90)
# # fig5_axes[0].set_title('Energy Level count given number of agents')
# # fig5_axes[0].set_xlabel('Site configuration')
# # fig5_axes[0].set_ylabel('Energy Level count')
# # # add text to the graph
# # for i, rect in enumerate(rects):
# #     height = rect.get_height()
# #     fig5_axes[0].text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),2)), ha='center', va='bottom')


# # possible_local_sites = calculate_possible_sites(global_params["influential_range"])
# # prob_dist_local_sites_given_motion = calculate_prob_dist_local_sites_given_motion(empty_motion_data, global_params);
# # possible_local_sites_data = np.array(list(zip(possible_local_sites,prob_dist_local_sites_given_motion)))
# # x_data = []
# # y_data = []
# # for local_sites_data in possible_local_sites_data:
# #     local_sites = local_sites_data[0]
# #     x_data.append(str(local_sites))
# #     y_data.append(calculate_entropy_for_local_sites_given_agent_motion(empty_motion_data, global_params, local_sites_data))
# # series_data = zip(x_data, y_data)
# # series_list = list(series_data)
# # series_list.sort(key=lambda x:x[1])
# # unzipped_data = list(zip(*series_list))
# # rects = fig5_axes[1].bar(unzipped_data[0], unzipped_data[1])
# # fig5_axes[1].tick_params('x', labelrotation=90)
# # fig5_axes[1].set_title('entropy of local sites given no recent motion')
# # fig5_axes[1].set_xlabel('Local Site configuration \n ')
# # fig5_axes[1].set_ylabel('Entropy')
# # # add text to the graph
# # for i, rect in enumerate(rects):
# #     height = rect.get_height()
# #     fig5_axes[1].text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),2)), ha='center', va='bottom')


# # fig5.tight_layout()







# # fig7, fig7_axes = plt.subplots(4, 1, constrained_layout=True)

# # possible_local_agents = [0,1,2,3]
# # x_data = []
# # y_data = []
# # # global_entropy = estimate_global_entropy(global_params)
# # for num_of_local_agents in possible_local_agents: 
# #     prob_dist = calculate_probs_dist_global_sites_given_local_agents(num_of_local_agents, global_params)
plt.savefig("prob_sytem_with_n_agents.svg")






#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################







# #     print(str(prob_dist))
#     # y_data.append(res)
#     # x_data.append(str(num_of_local_agents))
# # series_data = zip(x_data, y_data)
# # series_list = list(series_data)
# # series_list.sort(key=lambda x:x[1])
# # unzipped_data = list(zip(*series_list))
# # rects = fig7_axes.bar(unzipped_data[0], unzipped_data[1])
# # fig7_axes.tick_params('x', labelrotation=90)
# # fig7_axes.set_title('System entropy given number of local agents')
# # fig7_axes.set_xlabel('Number of agents in range')
# # fig7_axes.set_ylabel('Entropy')
# # # add text to the graph
# # for i, rect in enumerate(rects):
# #     height = rect.get_height()
# #     fig7_axes.text(rect.get_x() + rect.get_width()/2., 0.6*height,str(round(float(series_list[i][1]),2)), ha='center', va='bottom')

# # fig7.tight_layout()

# plt.show()








# fig = make_subplots(rows=3, cols=1)

# fig.append_trace(go.Scatter(
#     x=[3, 4, 5],
#     y=[1000, 1100, 1200],
# ), row=1, col=1)

# fig.append_trace(go.Scatter(
#     x=[2, 3, 4],
#     y=[100, 110, 120],
# ), row=2, col=1)

# fig.append_trace(go.Scatter(
#     x=[0, 1, 2],
#     y=[10, 11, 12]
# ), row=3, col=1)


# fig.update_layout(height=600, width=600, title_text="Stacked Subplots")
# fig.show()
