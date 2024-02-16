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

numOfAgentSet = [8, 16, 24, 32, 40, 48, 56, 64, 72];
interactionCoeff = 4.0;
externalFieldCoeff = 1.5;
lenOfTime = 64;

% Parameters for information quantification

influentialRange = 4;
timeWindowLength = 4;


%% This function calculates the probability of the local states through Monte Carlo Simulations

CAT = data_process_util;

for j = 1:9
metric_average_accumulate = 0.0;

numOfAgent = numOfAgentSet(j) 
    
for i = 1:1000
[Probability_Y_Given_SigmaM, interactionCoefficientVector, localStatesMatrix] = CAT.StatePredictionFromConditionalProbability(influentialRange, interactionCoeff, externalFieldCoeff);
% Calculate the time-space diagram
[Config, dataOfSpatialTemporal, StartPose] = CAT.StatisticalMechanicsBasedTraffic(numOfSite, numOfAgent, interactionCoefficientVector, externalFieldCoeff, lenOfTime, influentialRange);

% test of motion data
AgentIndex = 2;

t0 = 1;
seriesLenth = lenOfTime-timeWindowLength;
observabilityMetricSeries = zeros(seriesLenth, 1);
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
end

metric_average_accumulate = max(metric_average_accumulate, max(observabilityMetricSeries));
end
% test of motion data
metric_average_accumulate

end

x = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
y = [3, 4, 5, 6, 7];
z = [0.40573 0.51937 0.5558 0.61947 0.69629 0.74809 0.78989 0.80896 0.76330;
     0.5621 0.6472 0.7084 0.7667 0.8129 0.8513 0.8795 0.8888 0.8589;
     0.6602 0.7363 0.7995 0.8467 0.8843 0.9104 0.9268 0.9286 0.9066;
     0.7237 0.8023 0.8557 0.9005 0.9276 0.9436 0.9511 0.9495 0.9341;
     0.7731 0.8446 0.8968 0.9321 0.9513 0.9591 0.9632 0.9627 0.9521];

[M,c] = contour(x, y, z, 12, 'ShowText','on');
c.LineWidth = 3;
xlabel("Traffic Density");
ylabel("Influential Range");
zlabel("Entropy Based Observability");

