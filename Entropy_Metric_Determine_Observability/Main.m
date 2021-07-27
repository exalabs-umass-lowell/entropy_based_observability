% Readme:

% This file inclues an example of 8-site traffic system in a ring shape
% road, and analysis on its mutual information was conducted
% =========================================================================
%%  main function
clc;
clear;

% Traffic Model prameters
numOfSite = 64;
numOfAgent = 16;
interactionCoeff = 4.0;
externalFieldCoeff = 1.5;
lenOfTime = 64;

% Parameters for information quantification
influentialRange = 4;
timeWindowLength = 8;


%% This function calculates the probability of the local states through Monte Carlo Simulations

CAT = CATrafficDataProcess;
[Probability_Y_Given_SigmaM, interactionCoefficientVector, localStatesMatrix] = CAT.StatePredictionFromConditionalProbability(influentialRange, interactionCoeff, externalFieldCoeff);
% Calculate the time-space diagram
[Config, dataOfSpatialTemporal] = CAT.StatisticalMechanicsBasedTraffic(numOfSite, numOfAgent, interactionCoefficientVector, externalFieldCoeff, lenOfTime, influentialRange);
densityOfVehicle = numOfAgent/numOfSite;

imagesc(Config);
set(gca,'YDir','normal');
colormap(gray);
xlabel("Time (s/site)");
ylabel("Space (5 m/site)"); 

hold on;

% test of motion data
plot(1:lenOfTime, rem(cumsum(dataOfSpatialTemporal(2,:)),numOfSite+1), 'r', 'Linewidth',3);
xlabel("Time (1 s/step)");
ylabel("Space (5 m/site)");

t0 = 1;
historicalMotionData = dataOfSpatialTemporal(1, t0:t0+timeWindowLength-1);
ProbabilityOfFreeFlowFromData = CAT.MotionDataFiltering( historicalMotionData', timeWindowLength);
ProbabilityOfJamsFromData = 1-ProbabilityOfFreeFlowFromData;
% Probability_Y_Given_SigmaM  = CAT.StatePredictionFromConditionalProbability(influentialRange, interactionCoeff);
ProbabilityOfFreeFlowFromLocalConfig  = CAT.StatePredictionFromNeighboringSiteStates(influentialRange, Probability_Y_Given_SigmaM);
Probability_SigmaM_Given_Y  = CAT.DeriveAndNormalizeLikelihoodFromBayesian(Probability_Y_Given_SigmaM, ProbabilityOfFreeFlowFromLocalConfig, ProbabilityOfFreeFlowFromData);
[mutualInformation, observabilityMetric]  = CAT. observabilityQuantification(localStatesMatrix, Probability_SigmaM_Given_Y, numOfSite, numOfAgent, ProbabilityOfFreeFlowFromData);

%disp("Probability of the agent in free flow state: ");
%disp(ProbabilityOfFreeFlowFromData);

%disp("Probability of measurement given scenario: ");
%disp(Probability_Y_Given_SigmaM);

%disp("Probability of scenario given measurements: ");
%disp(Probability_SigmaM_Given_Y);

disp("Value of Mutual Information is calculated as: ");
disp(mutualInformation);

disp("Final results of Observability Metric is: ");
disp(observabilityMetric);

