% Readme:

% This file inclues an example of N-site traffic system in a ring shape
% road, and analysis focused on the observability variation at multiple
% spatial scales
% =========================================================================
%%  main function
clc;
clear;

% Traffic Model prameters
numOfSite = 64;
numOfAgent = 48;
interactionCoeff = 1.6;
externalFieldCoeff = 1.65;
lenOfTime = 64;
scalingConst = 2;


% Parameters for information quantification
influentialRange = 4;
timeWindowLength = 8;


% Scaling on arguments for traffic flow

numOfSite_3  = numOfSite/scalingConst^2;
numOfSite_2  = numOfSite/scalingConst;
numOfAgent_3 = numOfAgent/4;
initConfig_3 = [ones(numOfAgent_3,1);zeros(numOfSite_3-numOfAgent_3,1)];
initConfig_3 = initConfig_3(randperm(length(initConfig_3)));
initConfig_3 = initConfig_3*2-ones(numOfSite_3,1);  % changes states to be +/-1
initConfig_3(numOfSite_3*7/8:end, 1) = -1;
numOfAgent_3 = sum(initConfig_3 == 1);
numOfAgent_2 = 2*numOfAgent_3;
numOfAgent  =  4*numOfAgent_3;

initConfig_2 = reshape([initConfig_3'; initConfig_3'], [2*numOfSite_3, 1]);
initConfig   = reshape([initConfig_2'; initConfig_2'], [4*numOfSite_3, 1]);
lenOfTime_2 = lenOfTime/2;
lenOfTime_3 = lenOfTime_2/2;

influentialRange_2 = influentialRange/2;
influentialRange_3 = influentialRange_2/2;


CTM = CA_traffic_multiscale_modeling_util;
[interactionCoeff_2, externalFieldCoeff_2] = CTM.paramsRenormalization(interactionCoeff, externalFieldCoeff);
[interactionCoeff_3, externalFieldCoeff_3] = CTM.paramsRenormalization(interactionCoeff_2, externalFieldCoeff_2);
interactionCoefficientVector = CTM.HamiltonianInsideTrafficDynamics(influentialRange, interactionCoeff);% influentialRange, interactionCoeff, externalFieldCoeff);
interactionCoefficientVector_2 = CTM.HamiltonianInsideTrafficDynamics(influentialRange_2, interactionCoeff_2);% influentialRange, interactionCoeff, externalFieldCoeff);
interactionCoefficientVector_3 = CTM.HamiltonianInsideTrafficDynamics(influentialRange_3, interactionCoeff_3);% influentialRange, interactionCoeff, externalFieldCoeff);


[Config_ST, dataOfSpatialTemporal] = CTM.StatisticalMechanicsBasedTraffic(initConfig, interactionCoefficientVector, externalFieldCoeff, lenOfTime, influentialRange);

[Config_ST_2, dataOfSpatialTemporal_2] = CTM.StatisticalMechanicsBasedTraffic(initConfig_2, interactionCoefficientVector_2, externalFieldCoeff, lenOfTime_2, influentialRange_2);

[Config_ST_3, dataOfSpatialTemporal_3] = CTM.StatisticalMechanicsBasedTraffic(initConfig_3, interactionCoefficientVector_3, externalFieldCoeff, lenOfTime_3, influentialRange_3);






AgentIndex = 8;
AgentIndex_2 = AgentIndex/2;
AgentIndex_3 = AgentIndex/4;
%% Output the values/resutls

CAT = CA_traffice_data_process_util;
[Probability_Y_Given_SigmaM, interactionCoefficientVector, localStatesMatrix] = CAT.StatePredictionFromConditionalProbability(influentialRange, interactionCoeff, externalFieldCoeff);
[Config, dataOfSpatialTemp, StartPose] = CAT.StatisticalMechanicsBasedTraffic(numOfSite, numOfAgent, interactionCoefficientVector, externalFieldCoeff, lenOfTime, influentialRange);



t0 = 1;
seriesLenth = lenOfTime-timeWindowLength;
observabilityMetricSeries = zeros(seriesLenth, 1);
for t = t0:t0+seriesLenth-1
    historicalMotionData = dataOfSpatialTemporal(AgentIndex, t:t+timeWindowLength-1);
    ProbabilityOfFreeFlowFromData = CAT.MotionDataFiltering( historicalMotionData', timeWindowLength);
    ProbabilityOfJamsFromData = 1-ProbabilityOfFreeFlowFromData;
    % Probability_Y_Given_SigmaM  = CAT.StatePredictionFromConditionalProbability(influentialRange, interactionCoeff);
    [ProbabilityOfFreeFlowFromLocalConfig, ProbabilityOfLocalConfig]  = CAT.StatePredictionFromNeighboringSiteStates(influentialRange, Probability_Y_Given_SigmaM);
    Probability_SigmaM_Given_Y  = CAT.DeriveAndNormalizeLikelihoodFromBayesian(Probability_Y_Given_SigmaM, ProbabilityOfLocalConfig, ProbabilityOfFreeFlowFromData);
    entropyOfSigmaGlobalAtOriginalScale = CAT.SigmaGlobalCalculation(localStatesMatrix, numOfSite, numOfAgent);
    [mutualInformation, observabilityMetric]  = CAT.observabilityQuantification(localStatesMatrix, Probability_SigmaM_Given_Y, numOfSite, numOfAgent, ProbabilityOfFreeFlowFromData,  entropyOfSigmaGlobalAtOriginalScale);
    observabilityMetricSeries(t-t0+1, 1) = observabilityMetric;
       
end

%%


[Probability_Y_Given_SigmaM_2, interactionCoefficientVector, localStatesMatrix_2] = CAT.StatePredictionFromConditionalProbability(influentialRange_2, interactionCoeff_2, externalFieldCoeff_2);
[Config_2, dataOfSpatialTemp_2, StartPose_2] = CAT.StatisticalMechanicsBasedTraffic(numOfSite, numOfAgent, interactionCoefficientVector_2, externalFieldCoeff_2, lenOfTime_2, influentialRange_2);


t0 = 1;
observabilityMetricSeries_2 = zeros(seriesLenth, 1);

dataOfSpatialTemporal_2_new = zeros(1, size(dataOfSpatialTemporal_2,2)*2);
for i=2:size(dataOfSpatialTemporal_2, 2)
    if(dataOfSpatialTemporal_2(AgentIndex_2,i-1) == 0 && dataOfSpatialTemporal_2(AgentIndex_2,i) == 1)
        dataOfSpatialTemporal_2_new(1,2*i) = 1;
        if(randi(2) == 2)
            dataOfSpatialTemporal_2_new(1,2*i-1) = 0;
            dataOfSpatialTemporal_2_new(1,2*i-2)   = 1;
        else
            dataOfSpatialTemporal_2_new(1,2*i-2) = 1;
            dataOfSpatialTemporal_2_new(1,2*i-1)   = 0;
        end
    elseif(dataOfSpatialTemporal_2(AgentIndex_2,i-1) == 1 && dataOfSpatialTemporal_2(AgentIndex_2,i) == 1)
            dataOfSpatialTemporal_2_new(1,2*i-1) = 1;
            dataOfSpatialTemporal_2_new(1,2*i)   = 1;
    elseif(dataOfSpatialTemporal_2(AgentIndex_2,i) == 0)
            dataOfSpatialTemporal_2_new(1,2*i-1) = 0;
            dataOfSpatialTemporal_2_new(1,2*i)   = 0;
    end
end
pos_sequence_2 = cumsum(dataOfSpatialTemporal_2_new)+StartPose_2(AgentIndex_2);

for t = t0:t0+seriesLenth-1
    historicalMotionData = dataOfSpatialTemporal_2_new(1, t:t+timeWindowLength-1);
    ProbabilityOfFreeFlowFromData = CAT.MotionDataFiltering( historicalMotionData', timeWindowLength);
    ProbabilityOfJamsFromData = 1-ProbabilityOfFreeFlowFromData;
    % Probability_Y_Given_SigmaM  = CAT.StatePredictionFromConditionalProbability(influentialRange, interactionCoeff);
    [ProbabilityOfFreeFlowFromLocalConfig, ProbabilityOfLocalConfig]  = CAT.StatePredictionFromNeighboringSiteStates(influentialRange_2, Probability_Y_Given_SigmaM_2);
    Probability_SigmaM_Given_Y_2  = CAT.DeriveAndNormalizeLikelihoodFromBayesian(Probability_Y_Given_SigmaM_2, ProbabilityOfLocalConfig, ProbabilityOfFreeFlowFromData);
    [mutualInformation, observabilityMetric]  = CAT.observabilityQuantification(localStatesMatrix_2, Probability_SigmaM_Given_Y_2, numOfSite_2, numOfAgent_2, ProbabilityOfFreeFlowFromData, entropyOfSigmaGlobalAtOriginalScale);
    observabilityMetricSeries_2(t-t0+1, 1) = observabilityMetric;
       
end



%%
[Probability_Y_Given_SigmaM_3, interactionCoefficientVector, localStatesMatrix_3] = CAT.StatePredictionFromConditionalProbability(influentialRange_3, interactionCoeff_3, externalFieldCoeff_3);
[Config_3, dataOfSpatialTemp_3, StartPose_3] = CAT.StatisticalMechanicsBasedTraffic(numOfSite, numOfAgent, interactionCoefficientVector_3, externalFieldCoeff_3, lenOfTime_3, influentialRange_3);

dataOfSpatialTemporal_3_new = zeros(1, size(dataOfSpatialTemporal_3,2)*4);

for i=2:size(dataOfSpatialTemporal_3, 2)
    if(dataOfSpatialTemporal_3(AgentIndex_3,i-1) == 0 && dataOfSpatialTemporal_3(AgentIndex_3,i) == 1)
        dataOfSpatialTemporal_3_new(1, 4*i) = dataOfSpatialTemporal_3(AgentIndex_3,i);
        local_set = [1, 1, 0];
        dataOfSpatialTemporal_3_new(1,4*i-3:4*i-1) = local_set(randperm(3));
        dataOfSpatialTemporal_3_new(1,4*i) = 1;
    elseif(dataOfSpatialTemporal_3(AgentIndex_3,i-1) == 1 && dataOfSpatialTemporal_3(AgentIndex_3,i) == 1)
        dataOfSpatialTemporal_3_new(1,4*i-3:4*i) = [1,1,1,1];
    elseif( dataOfSpatialTemporal_3(AgentIndex_3,i-1) == 0 && dataOfSpatialTemporal_3(AgentIndex_3,i) == 0)
        dataOfSpatialTemporal_3_new(1,4*i-3:4*i) = [0, 0, 0, 0];
    else
        local_set = [1, 0, 0, 0];
        dataOfSpatialTemporal_3_new(1,4*i-3:4*i) = local_set(randperm(4));
    end
end
pos_sequence_3 = cumsum(dataOfSpatialTemporal_3_new)+StartPose_3(AgentIndex_3);


t0 = 1;

observabilityMetricSeries_3 = zeros(seriesLenth, 1);
for t = t0:t0+seriesLenth-1
    historicalMotionData = dataOfSpatialTemporal_3_new(1, t:t+timeWindowLength-1);
    ProbabilityOfFreeFlowFromData = CAT.MotionDataFiltering( historicalMotionData', timeWindowLength);
    % Probability_Y_Given_SigmaM  = CAT.StatePredictionFromConditionalProbability(influentialRange, interactionCoeff);
    [ProbabilityOfFreeFlowFromLocalConfig, ProbabilityOfLocalConfig]  = CAT.StatePredictionFromNeighboringSiteStates(influentialRange_3, Probability_Y_Given_SigmaM_3);
    Probability_SigmaM_Given_Y_3  = CAT.DeriveAndNormalizeLikelihoodFromBayesian(Probability_Y_Given_SigmaM_3, ProbabilityOfLocalConfig, ProbabilityOfFreeFlowFromData);
    [mutualInformation, observabilityMetric]  = CAT.observabilityQuantification(localStatesMatrix_3, Probability_SigmaM_Given_Y_3, numOfSite_3, numOfAgent_3, ProbabilityOfFreeFlowFromData,  entropyOfSigmaGlobalAtOriginalScale);
    observabilityMetricSeries_3(t-t0+1, 1) = observabilityMetric;
       
end


%%  Results visualization
figure
subplot(3,1,1)
imagesc(Config_ST);
set(gca,'YDir','normal');
colormap(gray);
hold on
% test of motion data
AgentIndex_ST = 8;
plot(1:lenOfTime, rem(cumsum(dataOfSpatialTemporal(AgentIndex_ST,:))+StartPose(AgentIndex_ST),numOfSite+1), 'r', 'LineWidth', 3);
xlabel("Time (s/site)");
ylabel("Space (5 m/site)"); 

subplot(3,1,2)
imagesc(Config_ST_2);
set(gca,'YDir','normal');
colormap(gray);
hold on
% test of motion data
AgentIndex_ST = 4;
plot(1:lenOfTime_2, rem(cumsum(dataOfSpatialTemporal_2(AgentIndex_ST,:))+StartPose_2(AgentIndex_ST),numOfSite+1), 'r', 'LineWidth', 3);
xlabel("Time (2 s/site)");
ylabel("Space (10 m/site)"); 

subplot(3,1,3)
imagesc(Config_ST_3);
set(gca,'YDir','normal');
colormap(gray);
hold on
% test of motion data
AgentIndex_ST = 2;
plot(1:lenOfTime_3, rem(cumsum(dataOfSpatialTemporal_3(AgentIndex_ST,:))+StartPose_3(AgentIndex_ST),numOfSite+1), 'r', 'LineWidth', 3);
xlabel("Time (4 s/site)");
ylabel("Space (20 m/site)"); 


%%
figure

subplot(3,1,1)
% test of motion data
yyaxis right
plot(1:lenOfTime, rem(cumsum(dataOfSpatialTemporal(AgentIndex,:))+StartPose(AgentIndex),numOfSite+1), 'r', 'LineWidth', 2);
ylabel("Space (5 m/site)");
yyaxis left
plot(timeWindowLength:seriesLenth+timeWindowLength-1, observabilityMetricSeries, 'k', 'LineWidth', 1);

xlabel("Time (s/timestep)");
ylabel("Observability Metric");
xlim([0 lenOfTime]);
ylim([0 1]);
grid on


subplot(3,1,2)
yyaxis right
plot(1:lenOfTime_2*2, pos_sequence_2, 'r', 'LineWidth', 2);
ylabel("Space (5 m/site)");
yyaxis left
%pos_sequence_2(isnan(pos_sequence_2))=1;
plot(timeWindowLength:seriesLenth+timeWindowLength-1, observabilityMetricSeries_2, 'k', 'LineWidth', 1);

xlabel("Time (s/timestep)");
ylabel("Observability Metric");
xlim([0 lenOfTime]);
ylim([0 1]);
grid on


subplot(3,1,3)
yyaxis right
plot(1:lenOfTime_3*4, pos_sequence_3, 'r', 'LineWidth', 2);
ylabel("Space (5 m/site)");
yyaxis left
%pos_sequence_3(isnan(pos_sequence_3))=1;
plot(timeWindowLength:seriesLenth+timeWindowLength-1, observabilityMetricSeries_3, 'k', 'LineWidth', 1);

xlabel("Time (s/timestep)");
ylabel("Observability Metric");
xlim([0 lenOfTime]);
ylim([0 1]);
grid on


