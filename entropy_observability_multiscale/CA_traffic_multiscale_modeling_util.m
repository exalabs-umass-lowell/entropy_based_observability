% Readme:

% This file inclues an example of N-site traffic system in a ring shape
% road, and analysis on its dynamical modeling at multiple spatial scales
% =========================================================================
%%  main function


classdef CA_traffic_multiscale_modeling_util
    
    methods
        %% This function determines the probability of measurements
        
        function [interactionCoefficient_scaled, externalFieldCoeff_scaled] = paramsRenormalization(obj, interactionCoefficient, externalFieldCoeff)
            
            interactionCoefficient_scaled = 1/2*log((exp(2*interactionCoefficient+2*externalFieldCoeff)+exp(-2*interactionCoefficient))/(exp(2*interactionCoefficient-2*externalFieldCoeff)+exp(-2*interactionCoefficient)));
            externalFieldCoeff_scaled = 1/4*log(((exp(2*interactionCoefficient+2*externalFieldCoeff)+exp(-2*interactionCoefficient))*(exp(2*interactionCoefficient-2*externalFieldCoeff)+exp(-2*interactionCoefficient)))/(exp(-2*externalFieldCoeff)+2+exp(2*externalFieldCoeff)));
                     
        end
        
               
        %% Interaction coefficient calculation
        
        function  interactionCoefficientVector = HamiltonianInsideTrafficDynamics(obj, influentialRange, interactionCoeff)% influentialRange, interactionCoeff, externalFieldCoeff)
            
            interactionCoefficientVector = zeros(influentialRange, 1);
            for index = 1:influentialRange
               interactionCoefficientVector(index) = tanh(interactionCoeff)^index;
            end
        end
        
        
        %% This function would determine the dynamics of CA traffic model
        
        function [Config, motionRecord] = StatisticalMechanicsBasedTraffic(obj, initConfig, interactionCoefficientVector, externalFieldCoeff, lenOfTime, influentialRange)

            numOfSite = size(initConfig,1);
            numOfAgent = sum(initConfig == 1);
            % Scaling over the sensor range and road length
            Config = -ones(numOfSite, lenOfTime);
            Config(:,1) = initConfig;
            motionRecord = zeros(numOfAgent, lenOfTime);
            localStates = zeros(influentialRange, 1);
            Hash_Pose2Agent = zeros(numOfSite, lenOfTime);
            
            AgentInd = 1;
            for Ind = 1:numOfSite
               if (initConfig(Ind)== 1)
                  Hash_Pose2Agent(Ind,1) = AgentInd;
                  AgentInd = AgentInd+1;
               end
            end
            
            for t = 1:lenOfTime-1
                i = 1;
                while(i<=numOfSite)
                    % In case get to the boundary
                    ind = i;
                    ind2 = i+1;
                    if ind == numOfSite
                        ind2 = 1;
                    end
                    
                    for iter = ind:ind+influentialRange-1
                       if(iter > numOfSite)
                           localStates(iter-ind+1) =  Config(iter-numOfSite+1,t);
                       else
                           localStates(iter-ind+1) = Config(iter,t);
                       end
                    end
                    % Interaction terms
                    localStates = (localStates==1).*1;
                    H_int = interactionCoefficientVector'*localStates;
                    % External field term
                    H_ext = externalFieldCoeff*Config(ind,t);
                    % Metroplis Algorithm to update the sites
                    H = H_int - H_ext;
                    P_tr = min(1, exp(-H));               % probability based on Hamiltonian
                    Rand_Val = rand;             % produce a random number
       
                    % make sure the vehicle are moving forward/backward
                    if Config(ind, t) == 1 && Config(ind2,t)==-1 && Rand_Val <= P_tr
                        Config(ind,t+1) = Config(ind2,t);
                        Config(ind2,t+1) = Config(ind,t);
                        i = i+1;
                        AgentIndex = Hash_Pose2Agent(ind,t);                        
                        motionRecord(AgentIndex, t+1) = 1;
                        Hash_Pose2Agent(ind2,t+1) = AgentIndex;
                    else
                        Config(ind,t+1) = Config(ind,t);
                        Hash_Pose2Agent(ind,t+1) = Hash_Pose2Agent(ind,t);
                    end
                    i = i+1;
                end
            end
        end
        
        % ==========================================================

    end
end