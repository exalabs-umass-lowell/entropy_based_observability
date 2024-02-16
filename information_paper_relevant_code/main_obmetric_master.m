% Readme:

% This file inclues an example of N-site traffic system in a ring shape
% road, and analysis on its mutual information was conducted to extract
% observability metric
% =========================================================================
%%  main function
clc;
clear;

% Traffic Model prameters
numOfSite = 80;
numOfAgent = 8;
interactionCoeff = 4.0;
externalFieldCoeff = 1.5;
lenOfTime = 64;

% Parameters for information quantification
influentialRange = 3;
timeWindowLength = 4;


%% This function calculates the probability of the local states through Monte Carlo Simulations

CAT = data_process_util;
[Probability_Y_Given_SigmaM, interactionCoefficientVector, localStatesMatrix] = CAT.StatePredictionFromConditionalProbability(influentialRange, interactionCoeff, externalFieldCoeff);
% Calculate the time-space diagram
[Config, dataOfSpatialTemporal, StartPose] = CAT.StatisticalMechanicsBasedTraffic(numOfSite, numOfAgent, interactionCoefficientVector, externalFieldCoeff, lenOfTime, influentialRange);
densityOfVehicle = numOfAgent/numOfSite;


imagesc(Config);
set(gca,'YDir','normal');
colormap(gray);
xlabel("Time (s/site)");
ylabel("Space (5 m/site)");

hold on;

% test of motion data
AgentIndex = 3;
plot(1:lenOfTime, rem(cumsum(dataOfSpatialTemporal(AgentIndex,:))+StartPose(AgentIndex),numOfSite+1), 'r', 'LineWidth', 3);
xlabel("Time (s/timestep)");
ylabel("Space (5 m/site)");


t0 = 1;
seriesLenth = lenOfTime-timeWindowLength;
observabilityMetricSeries = zeros(seriesLenth, 1);

% motion_data = [0, 0, 0, 0; 0, 0, 0, 1; 0, 0, 1, 0; 0, 0, 1, 1; 0, 1, 0, 0;
%                0, 1, 0, 1; 0, 1, 1, 0; 0, 1, 1, 1; 1, 0, 0, 0; 1, 0, 0, 1;
%                1, 0, 1, 0; 1, 0, 1, 1;1, 1, 0, 0; 1, 1, 0, 1; 1, 1, 1, 0;
%                1, 1, 1, 1];

entropy_data = zeros(16, 1);

for t = t0:t0+seriesLenth-1
    historicalMotionData = dataOfSpatialTemporal(AgentIndex, t:t+timeWindowLength-1);
    % historicalMotionData = [ones(1, timeWindowLength-2), 1, 1];
    % historicalMotionData = historicalMotionData(index,:);
    ProbabilityOfFreeFlowFromData = CAT.MotionDataFiltering(historicalMotionData', timeWindowLength);
    ProbabilityOfJamsFromData = 1-ProbabilityOfFreeFlowFromData;
    Probability_Y_Given_SigmaM  = CAT.StatePredictionFromConditionalProbability(influentialRange, interactionCoeff, externalFieldCoeff);
    [ProbabilityOfFreeFlowFromLocalConfig, ProbabilityOfLocalConfig]  = CAT.StatePredictionFromNeighboringSiteStates(influentialRange, Probability_Y_Given_SigmaM);
    Probability_SigmaM_Given_Y  = CAT.DeriveAndNormalizeLikelihoodFromBayesian(Probability_Y_Given_SigmaM, ProbabilityOfLocalConfig, ProbabilityOfFreeFlowFromData);
    entropyOfSigmaGlobalAtOriginalScale = CAT.SigmaGlobalCalculation(localStatesMatrix, numOfSite, numOfAgent);
    
    [mutualInformation, observabilityMetric, sumOfConditionalEntropy_SigmaGlobal_Given_Y]  = CAT.observabilityQuantification(localStatesMatrix, Probability_SigmaM_Given_Y, numOfSite, numOfAgent, ProbabilityOfFreeFlowFromData, entropyOfSigmaGlobalAtOriginalScale);
    % entropy_data(index) = sumOfConditionalEntropy_SigmaGlobal_Given_Y;
    observabilityMetricSeries(t-t0+1, 1) = observabilityMetric;
    
    
    %disp("Value of Mutual Information is calculated as: ");
    %disp(mutualInformation);
    
    %disp("Final results of Observability Metric is: ");
    %disp(observabilityMetric);    
end

figure
% test of motion data
yyaxis right
plot(1:lenOfTime, rem(cumsum(dataOfSpatialTemporal(AgentIndex,:))+StartPose(AgentIndex),numOfSite+1), 'r', 'LineWidth', 2);
ylabel("Space (5 m/site)");
yyaxis left
plot(timeWindowLength:seriesLenth+timeWindowLength-1, observabilityMetricSeries, 'k', 'LineWidth', 1);

xlabel("Time (s/timestep)");
ylabel("Observability Metric");
ylim([0 1]);
grid on


% figure
% test of motion data
% yyaxis right
% plot(1:lenOfTime, rem(cumsum(dataOfSpatialTemporal(AgentIndex,:))+StartPose(AgentIndex),numOfSite+1), 'r', 'LineWidth', 2);
% ylabel("Space (5 m/site)");
% yyaxis left
% plot(timeWindowLength:seriesLenth+timeWindowLength-2, observabilityMetricSeries(2:end)-observabilityMetricSeries(1:end-1), 'k', 'LineWidth', 1);
% 
% xlabel("Time (s/timestep)");
% ylabel("Observability Metric Variation Rate");
% ylim([-0.2 0.2]);
% grid on


