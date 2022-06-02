# system libs
import sys
import math

# dependencies - numpy, matplotlib
import numpy as np 
from numpy.random import rand
from numpy import log, dot, e
import matplotlib.pyplot as plt

# local files
from utility import *
from simulation import *
from analysis import *

np.set_printoptions(threshold=sys.maxsize)

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

# calculate the time-space diagram
simulation_data = statistical_mechanics_based_traffic(global_params)

# extract simulation data
sites_through_time = simulation_data["sites"]
sites_agent_occupancy_through_time = simulation_data["sites_agent_occupancy"]
agents_motion_through_time = simulation_data["agents_motion"]

agents_start_positions = get_agent_start_positions(sites_agent_occupancy_through_time, global_params)
densityOfVehicle = global_params["num_of_agents"]/global_params["num_of_sites"];

# test of motion data
AgentIndex = 0;
# plt.plot(np.arange(0, global_params["len_of_time"]), np.remainder(np.cumsum(dataOfSpatialTemporal[AgentIndex,:])+StartPose[AgentIndex],global_params["num_of_sites"]+1), 'r', 'LineWidth', 3)
# plt.xlabel("Time (s/timestep)")
# plt.ylabel("Space (5 m/site)")

# test of motion data
Agent2Index = 1;
# plt.plot(np.arange(0, global_params["len_of_time"]), np.remainder(np.cumsum(dataOfSpatialTemporal[Agent2Index,:])+StartPose[Agent2Index],global_params["num_of_sites"]+1), 'r', 'LineWidth', 3);
# plt.xlabel("Time (s/timestep)");
# plt.ylabel("Space (5 m/site)");

# prob_y_given_sigma_M, local_states_mat = state_prediction_from_conditional_prob(global_params);

seriesLenth = global_params["len_of_time"]-global_params["len_of_time_window"]
observability_series = get_observability_metric_series_for_agent(agents_motion_through_time, global_params, AgentIndex)
observability_series_2 = get_observability_metric_series_for_agent(agents_motion_through_time, global_params, Agent2Index)

# figure
fig, ax_left = plt.subplots(1,1)
ax_right = ax_left.twinx()
# ax_left.plot(np.arange(0,global_params["len_of_time"]), np.remainder(np.cumsum(agents_motion_through_time[AgentIndex,:])+agents_start_positions[AgentIndex],global_params["num_of_sites"]+1), '--r', 'LineWidth', 3)
x_series = np.remainder(np.cumsum(agents_motion_through_time[AgentIndex,:])+agents_start_positions[AgentIndex],global_params["num_of_sites"]+1)
ax_left.plot(np.arange(0,global_params["len_of_time"]), x_series, '--r', 'LineWidth', 3)
ax_left.plot(np.arange(0,global_params["len_of_time"]), np.remainder(np.cumsum(agents_motion_through_time[Agent2Index,:])+agents_start_positions[Agent2Index],global_params["num_of_sites"]+1), '--g', 'LineWidth', 3)
# ax_left.plot(np.arange(0,global_params["len_of_time"]), np.remainder(np.cumsum(agents_motion_through_time[Agent2Index,:])+agents_start_positions[Agent2Index],global_params["num_of_sites"]+1), '--g', 'LineWidth', 3)
ax_left.set_ylabel("Space (5 m/site)");
# these show the observability starting at the time window length then going to the end (if the time window measurement is taken from the past)
ax_right.plot(np.arange(global_params["len_of_time_window"], seriesLenth+global_params["len_of_time_window"]), observability_series, '-or', 'LineWidth', 3)
ax_right.plot(np.arange(global_params["len_of_time_window"], seriesLenth+global_params["len_of_time_window"]), observability_series_2, '-og', 'LineWidth', 3)
# these show the observability starting at the begining then going to one time window length befor the end (if the time window measuerment is taken from the future)
# ax_right.plot(np.arange(0, global_params["len_of_time"]- global_params["len_of_time_window"]), observability_series, '-or', 'LineWidth', 3)
# ax_right.plot(np.arange(0, global_params["len_of_time"]- global_params["len_of_time_window"]), observability_series_2, '-og', 'LineWidth', 3)
ax_right.set_ybound(0,1)

plt.xlabel("Time (s/timestep)");
ax_right.set_ylabel("Observability Metric");
plt.show();
   