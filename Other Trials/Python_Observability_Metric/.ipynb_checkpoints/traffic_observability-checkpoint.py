import numpy as np 
from numpy.random import rand
from numpy import log, dot, e

import matplotlib.pyplot as plt


# % Traffic Model prameters
numOfSite = 64;
numOfAgent = 16;
interaction_coeff = 4.0;
external_field_coeff = 1.5;
lenOfTime = 64;

# % Parameters for information quantification
influential_range = 4;
timeWindowLength = 16;


def de2bi(d):
    n = len(d)
    power = 2**np.arange(n)
    d = d * np.ones((1,n))
    b = np.floor((d%(2*power))/power)
    return b


# %% This function helps determine the probability of moving forward through the neighboring site states within influential range     
# function [Probability_Y_Given_SigmaM, interactionCoefficientVector, localStatesMatrix] = StatePredictionFromConditionalProbability(obj, influentialRange, InteractionCoeff, externalFieldCoeff)
def state_prediction_from_conditional_prob(influential_range, interaction_coeff, external_field_coeff):
    num_of_local_config = 2^influential_range;
    #Probability_Y_Given_SigmaM = zeros(num_of_local_config, 1);
    prob_y_given_sigma_M = np.zeros((num_of_local_config, 1));

    # tempStateMatrix = de2bi(1:num_of_local_config);
    state_range = np.arange(1, num_of_local_config);
    temp_state_mat = np.zeros((state_range, 1));

    # local_states_mat = tempStateMatrix(:, 1:influential_range);
    influential_range_vec = np.arange(1, influential_range)
    local_states_mat = temp_state_mat[:, influential_range_vec];

    # interaction_coeff_vec = zeros(influential_range, 1);
    interaction_coeff_vec = np.zeros((influential_range, 1));
    # for index = 1:influential_range:
    for i in influential_range_vec:
        interaction_coeff_vec[i] = np.tanh(interaction_coeff)**i;

    # % localStates = (localStates==1).*1;    % convert to 1/0 vector
    # for index = 1:num_of_local_config:
    for i in state_range:
        #Probability_Y_Given_SigmaM(index) = min(1, exp(external_field_coeff - interaction_coeff_vec'*local_states_mat(index,:)'));
        prob_y_given_sigma_M[i] = min(1, np.exp(external_field_coeff - interaction_coeff_vec.conjugate().transpose() * local_states_mat[i,:].conjugate().transpose()));

    return prob_y_given_sigma_M, interaction_coeff_vec, local_states_mat


prob_y_given_sigma_M, interaction_coeff_vec, local_states_mat = state_prediction_from_conditional_prob(influential_range, interaction_coeff, external_field_coeff);


# function [Config, motionRecord, startPose] = StatisticalMechanicsBasedTraffic(obj, numOfSite, numOfAgent, interactionCoefficientVector, externalFieldCoeff, lenOfTime, influentialRange)
            
#             initConfig = [ones(numOfAgent,1);zeros(numOfSite-numOfAgent,1)];
#             initConfig = initConfig(randperm(length(initConfig)));
#             initConfig = initConfig*2-ones(numOfSite,1);  % changes states to be +/-1
            
#             % Scaling over the sensor range and road length
#             Config = -ones(numOfSite, lenOfTime);
#             Config(:,1) = initConfig;
#             motionRecord = zeros(numOfAgent, lenOfTime);
#             startPose = zeros(numOfAgent, 1);
#             localStates = zeros(influentialRange, 1);
#             Hash_Pose2Agent = zeros(numOfSite, lenOfTime);
            
#             AgentInd = 1;
#             for Ind = 1:numOfSite
#                 if (initConfig(Ind)== 1)
#                     Hash_Pose2Agent(Ind,1) = AgentInd;
#                     AgentInd = AgentInd+1;
#                 end
#             end
            
#             for t = 1:lenOfTime-1
#                 i = 1;
#                 while(i<=numOfSite)
#                     % In case get to the boundary
#                     ind = i;
#                     ind2 = i+1;
#                     if ind == numOfSite
#                         ind2 = 1;
#                     end
                    
#                     for iter = ind:ind+influentialRange-1
#                         if(iter > numOfSite)
#                             localStates(iter-ind+1) =  Config(iter-numOfSite+1,t);
#                         else
#                             localStates(iter-ind+1) = Config(iter,t);
#                         end
#                     end
#                     % Interaction terms
#                     localStates = (localStates==1).*1;
#                     H_int = interactionCoefficientVector'*localStates;
#                     % External field term
#                     H_ext = externalFieldCoeff*Config(ind,t);
#                     % Metroplis Algorithm to update the sites
#                     H = H_int - H_ext;
#                     P_tr = min(1, exp(-H));               % probability based on Hamiltonian
#                     Rand_Val = rand;             % produce a random number
                    
#                     if t == 1 && Config(ind, t) == 1
#                         AgentIndex = Hash_Pose2Agent(ind,t);
#                         startPose(AgentIndex,1) = ind;
#                     end
#                     % make sure the vehicle are moving forward/backward
#                     if Config(ind, t) == 1 && Config(ind2,t)==-1 && Rand_Val <= P_tr
#                         Config(ind,t+1) = Config(ind2,t);
#                         Config(ind2,t+1) = Config(ind,t);
#                         i = i+1;
#                         AgentIndex = Hash_Pose2Agent(ind,t);
#                         motionRecord(AgentIndex, t+1) = 1;
                        
#                         Hash_Pose2Agent(ind2,t+1) = AgentIndex;
#                     else
#                         Config(ind,t+1) = Config(ind,t);
#                         Hash_Pose2Agent(ind,t+1) = Hash_Pose2Agent(ind,t);
#                     end
#                     i = i+1;
#                 end
#             end
#         end

