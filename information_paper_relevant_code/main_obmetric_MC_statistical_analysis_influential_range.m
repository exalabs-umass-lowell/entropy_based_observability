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

influentialRange = 7;
timeWindowLength = 4;


%% This function calculates the probability of the local states through Monte Carlo Simulations

CAT = data_process_util;

for j = 1:9
metric_average_accumulate = 0.0;

numOfAgent = numOfAgentSet(j) 
    
for i = 1:100
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

metric_average_accumulate = metric_average_accumulate + mean(observabilityMetricSeries);
end
% test of motion data
metric_average_accumulate/100

end



