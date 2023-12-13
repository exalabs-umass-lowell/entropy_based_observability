% Readme:

% This file inclues an example of 8-site traffic system in a ring shape
% road, and analysis on its mutual information was conducted
% =========================================================================
%%  main function


classdef CATrafficDataProcess
    
    methods
        %% This function determines the probability at the states of being moving foward during measurements Y
        
        function ProbabilityOfFreeFlowFromData = MotionDataFiltering(obj, historicalMotionData, timeWindowLength)
            % create a gaussian filter to process the data
            gaussianSequence = gausswin(timeWindowLength*2);
            gaussianKernel = gaussianSequence(1:timeWindowLength)/sum(gaussianSequence)*2;
            ProbabilityOfFreeFlowFromData = gaussianKernel'*historicalMotionData;
        end
        
        % ==========================================================
        
        %% This function helps determine the probability of moving forward through the neighboring site states within influential range
             
        function [Probability_Y_Given_SigmaM, interactionCoefficientVector, localStatesMatrix]  = StatePredictionFromConditionalProbability(obj, influentialRange, InteractionCoeff, externalFieldCoeff)
            NumOfLocalConfig = 2^influentialRange;
            Probability_Y_Given_SigmaM = zeros(NumOfLocalConfig, 1);
            tempStateMatrix = de2bi(1:NumOfLocalConfig);
            localStatesMatrix = tempStateMatrix(:, 1:influentialRange);
            interactionCoefficientVector = zeros(influentialRange, 1);
            for index = 1:influentialRange
                interactionCoefficientVector(index) = tanh(InteractionCoeff)^index;
            end
            % localStates = (localStates==1).*1;    % convert to 1/0 vector
            for index = 1:NumOfLocalConfig
                Probability_Y_Given_SigmaM(index) = min(1, exp(-externalFieldCoeff -interactionCoefficientVector'*localStatesMatrix(index,:)'));
            end
        end
        
        % ==========================================================
        
        %% This function helps determine the probability of the local states
        
        function [ProbabilityOfFreeFlowFromLocalConfig, ProbabilityOfLocalConfig]  = StatePredictionFromNeighboringSiteStates(obj, influentialRange, Probability_Y_Given_SigmaM)
            NumOfLocalConfig = 2^influentialRange;
            ProbabilityOfLocalConfig = zeros(NumOfLocalConfig, 1);
            
            for index = 1:NumOfLocalConfig
                ProbabilityOfLocalConfig(index) = 1/NumOfLocalConfig;
            end
            ProbabilityOfFreeFlowFromLocalConfig = ProbabilityOfLocalConfig'*Probability_Y_Given_SigmaM;
        end
        
        % ==========================================================
        
        %% This function determines the likelihood of local states given measurement and conduct normalization
        
        function Probability_SigmaM_Given_Y  = DeriveAndNormalizeLikelihoodFromBayesian(obj, Probability_Y_Given_SigmaM, ProbabilityOfLocalConfig, ProbabilityOfFreeFlowFromData)
            % Calculation of likelihood
            if ProbabilityOfFreeFlowFromData == 1
                Likelihood_SigmaM_Given_Y   = Probability_Y_Given_SigmaM.*ProbabilityOfLocalConfig/ProbabilityOfFreeFlowFromData;
            elseif ProbabilityOfFreeFlowFromData == 0
                Likelihood_SigmaM_Given_Y  = (1 - Probability_Y_Given_SigmaM).*ProbabilityOfLocalConfig/(1-ProbabilityOfFreeFlowFromData);
            else
                Likelihood_SigmaM_Given_Y   = Probability_Y_Given_SigmaM.*ProbabilityOfLocalConfig/ProbabilityOfFreeFlowFromData + (1 - Probability_Y_Given_SigmaM).*ProbabilityOfLocalConfig/(1-ProbabilityOfFreeFlowFromData);
            end
            % Conduct normalization on likelihood to derive the probability
            Probability_SigmaM_Given_Y  = Likelihood_SigmaM_Given_Y/sum(Likelihood_SigmaM_Given_Y(:));
        end
        % ==========================================================
        
        

        
        %% This function calculates the entropy of global systems at original scale
        function entropyOfSigmaGlobal = SigmaGlobalCalculation(obj, localStates, numOfSite, numOfAgent)
            
            influentialRange = size(localStates, 2);
            entropyOfSigmaGlobal = 0;
            for agentNum = 1:min(numOfAgent, influentialRange)
            % Calculation of gloabl state entropy   H(\Simga)
                numOfAgentOutOfRange = numOfAgent - agentNum;
                numOfLocalConfig = nchoosek(influentialRange, agentNum);
                probabilityOfLocalConfig = 1/numOfLocalConfig;
                energyLevel = obj.energyLevelCalculation(numOfAgentOutOfRange, numOfSite-1-influentialRange);
                
                probabilityGlobalDistributionMatrix = energyLevel/sum(energyLevel(:))*probabilityOfLocalConfig*ones(numOfLocalConfig,1)';
                probabilityGlobalDistributionVec = reshape(probabilityGlobalDistributionMatrix',1,[]);
                entropyOfSigmaGlobal = entropyOfSigmaGlobal + obj.entropyCalculation(probabilityGlobalDistributionVec);  %H(\Sigma)
            end
                    
        end 
            
        %% This function helps determine the probability of the global states
        % localStates:   cover all possible config within influentia
        % range        #_config  *  influential_range
        
        function [mutualInformation, observabilityMetric, sumOfConditionalEntropy_SigmaGlobal_Given_Y] = observabilityQuantification(obj, localStates, Probability_SigmaM_Given_Y, numOfSite, numOfAgent, ProbabilityOfFreeFlowFromData, entropyOfSigmaGlobalAtOriginalScale)
            numOfAgentInfluenced = sum(localStates, 2);
            numOfLocalConfig = size(localStates, 1);
            influentialRange = size(localStates, 2);
            numOfAgentOutOfRange = numOfAgent - numOfAgentInfluenced;
            sumOfConditionalEntropy_SigmaGlobal_Given_Y = 0;
            
            
            % iterate all possible configuration
            % ==================================================================================================
            % Calculation of conditional entropy   H(\Sigma|Y)
            for index = 1:numOfLocalConfig
                %  Calculation the energy distributions of sites beyond
                %  influential range   P(\Sigma_complementary | \Sigma_m)
                energyLevel = obj.energyLevelCalculation(numOfAgentOutOfRange(index), numOfSite-1-influentialRange);
                Probability_SigmaComp_Given_SigmaM = energyLevel/sum(energyLevel(:));
                %  Calculate the probability distributions by merging
                %  distributions within and beyond influential range
                %  P(\Sigma_complementary | \Sigma_m)* P(\Sigma_m| Y)
                Probability_SigmaGlobal_Given_Y = Probability_SigmaComp_Given_SigmaM*Probability_SigmaM_Given_Y(index);
                sumOfConditionalEntropy_SigmaGlobal_Given_Y  = sumOfConditionalEntropy_SigmaGlobal_Given_Y + ProbabilityOfFreeFlowFromData*obj.entropyCalculation(Probability_SigmaGlobal_Given_Y) + ...
                    (1-ProbabilityOfFreeFlowFromData)*obj.entropyCalculation(1-Probability_SigmaGlobal_Given_Y);
            end
            %====================================================================================================
            
            % definition of mutual information

            mutualInformation = entropyOfSigmaGlobalAtOriginalScale - sumOfConditionalEntropy_SigmaGlobal_Given_Y;
            % probabilityMeasurementDistribution = [ProbabilityOfFreeFlowFromData, 1-ProbabilityOfFreeFlowFromData];
            % entropyOfMeasurement = obj.entropyCalculation(probabilityMeasurementDistribution);   %H(Y)
            entropyOfSigmaGlobal = obj.SigmaGlobalCalculation(localStates,  numOfSite, numOfAgent);
            observabilityMetric = mutualInformation/entropyOfSigmaGlobal;
            % measurementEfficiencyMetric = mutualInformation/entropyOfMeasurement;
        end
        
        % ============================================================
        
        %% This function helps determine the probability of the global states
        
        function  countOfEnergyLevel = energyLevelCalculation(obj, numOfAgent, numOfSite)
            if(numOfSite>=2*numOfAgent)
                energyLevel = (numOfSite-2*numOfAgent-1)-2*numOfAgent:4:numOfSite-1-2;
            else
                energyLevel = (2*numOfAgent-numOfSite-1)-2*(numOfSite-numOfAgent):4:numOfSite-1-2;
            end
            numOfLevel = size(energyLevel,2);
            numOfVacantSite = numOfSite - numOfAgent;
            countOfEnergyLevel = zeros(numOfLevel,1);
            %   find spaces between vacant sites and insert set of vehicles  *  group vehicles into (index-1) sets
            if(numOfVacantSite >= numOfAgent)
                for index = 1:numOfLevel
                    countOfEnergyLevel(index) = nchoosek(numOfVacantSite, numOfAgent+1-index)*nchoosek(numOfAgent-1, index-1);
                end
            else
                 for index = 1:numOfLevel
                    countOfEnergyLevel(index) = nchoosek(numOfAgent, numOfVacantSite+1-index)*nchoosek(numOfVacantSite-1, index-1);
                end
            end
        end
        % ============================================================
        
        
        %% This function helps determine the probability of the global states
        
        function res  = entropyCalculation(obj,probabilityDistribution)
            probabilityDistribution_m = nonzeros(probabilityDistribution);
            res = 0;
            totalLength = size(probabilityDistribution_m,1);
            for i = 1:totalLength
                res = res - probabilityDistribution_m(i)*log(probabilityDistribution_m(i));
            end
        end
        % ============================================================
        
        
        
        
        %% This function would determine the dynamics of CA traffic model
        
        function [Config, motionRecord, startPose, p_sequence] = StatisticalMechanicsBasedTraffic(obj, numOfSite, numOfAgent, interactionCoefficientVector, externalFieldCoeff, lenOfTime, influentialRange)
            
            initConfig = [ones(numOfAgent,1);zeros(numOfSite-numOfAgent,1)];
            initConfig = initConfig(randperm(length(initConfig)));
            initConfig = initConfig*2-ones(numOfSite,1);  % changes states to be +/-1
            
            % Scaling over the sensor range and road length
            Config = -ones(numOfSite, lenOfTime);
            Config(:,1) = initConfig;
            motionRecord = zeros(numOfAgent, lenOfTime);
            startPose = zeros(numOfAgent, 1);
            localStates = zeros(influentialRange, 1);
            Hash_Pose2Agent = zeros(numOfSite, lenOfTime);
            
            AgentInd = 1;
            for Ind = 1:numOfSite
                if (initConfig(Ind)== 1)
                    Hash_Pose2Agent(Ind,1) = AgentInd;
                    AgentInd = AgentInd+1;
                end
            end
            
            p_sequence = zeros(lenOfTime, 1);
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
                    p_sequence(t) = P_tr;
                    if t == 1 && Config(ind, t) == 1
                        AgentIndex = Hash_Pose2Agent(ind,t);
                        startPose(AgentIndex,1) = ind;
                    end
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