% Readme:

% This file inclues an example of N-site traffic system in a ring shape
% road, and analysis on its mutual information was conducted to extract
% observability metric
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
timeWindowSet = [8, 12, 16];
AgentIndex = 3;


%% This function calculates the probability of the local states through Monte Carlo Simulations

CAT = CATrafficDataProcess;
[Probability_Y_Given_SigmaM, interactionCoefficientVector, localStatesMatrix] = CAT.StatePredictionFromConditionalProbability(influentialRange, interactionCoeff, externalFieldCoeff);
% Calculate the time-space diagram
[Config, dataOfSpatialTemporal, StartPose] = CAT.StatisticalMechanicsBasedTraffic(numOfSite, numOfAgent, interactionCoefficientVector, externalFieldCoeff, lenOfTime, influentialRange);
densityOfVehicle = numOfAgent/numOfSite;



%%  Visualize on the variation of observability metric according to different time window
t0 = 1;
seriesLenth = lenOfTime-timeWindowSet(1);
observabilityMetricSequence = zeros(seriesLenth, size(timeWindowSet,2));
mutualInformationSequence = zeros(seriesLenth, size(timeWindowSet,2));

for index = 1:size(timeWindowSet,2)    
    timeWindowLength = timeWindowSet(index);
    for t = t0:t0+lenOfTime-timeWindowLength-1
        historicalMotionData = dataOfSpatialTemporal(AgentIndex, t:t+timeWindowLength-1);
        ProbabilityOfFreeFlowFromData = CAT.MotionDataFiltering( historicalMotionData', timeWindowLength);
        ProbabilityOfJamsFromData = 1-ProbabilityOfFreeFlowFromData;
        [ProbabilityOfFreeFlowFromLocalConfig, ProbabilityOfLocalConfig]  = CAT.StatePredictionFromNeighboringSiteStates(influentialRange, Probability_Y_Given_SigmaM);
        Probability_SigmaM_Given_Y  = CAT.DeriveAndNormalizeLikelihoodFromBayesian(Probability_Y_Given_SigmaM, ProbabilityOfLocalConfig, ProbabilityOfFreeFlowFromData);
        [mutualInformation, observabilityMetric]  = CAT.observabilityQuantification(localStatesMatrix, Probability_SigmaM_Given_Y, numOfSite, numOfAgent, ProbabilityOfFreeFlowFromData);
        observabilityMetricSequence(t-t0+1, index) = observabilityMetric;
        mutualInformationSequence(t-t0+1, index)= mutualInformation;
    end
end

figure
% test of motion data
yyaxis right
plot(1:lenOfTime, rem(cumsum(dataOfSpatialTemporal(AgentIndex,:))+StartPose(AgentIndex),numOfSite+1), 'r', 'LineWidth', 2);
ylabel("Space (5 m/site)");
yyaxis left
h = plot(timeWindowSet(1):seriesLenth+timeWindowSet(1)-1, observabilityMetricSequence(:,1), 'k', ... 
    timeWindowSet(2):seriesLenth+timeWindowSet(2)-1, observabilityMetricSequence(:,2), 'b', ...
    timeWindowSet(3):seriesLenth+timeWindowSet(3)-1, observabilityMetricSequence(:,3), 'g');
set(h, {'LineWidth'}, {2;2;2});
xlabel("Time (s/timestep)");
ylabel("Observability Metric");
legend('TimeWindow = 8', 'TimeWindow=12', 'TimeWindow=16');
xlim([0 62]);
ylim([0 1]);
grid on




%%  Visualize on the variation of entropyOfSigmaGlobalSequence && ConditionalEntropy_SigmaGlobal_Given_Y

entropyOfSigmaGlobalSequence = mutualInformationSequence(:,1)./observabilityMetricSequence(:,1);
conditionalEntropyOfSigmaGlobalGivenYSequence = entropyOfSigmaGlobalSequence - mutualInformationSequence(:,1);

figure
% test of motion data
yyaxis right
plot(1:lenOfTime, rem(cumsum(dataOfSpatialTemporal(AgentIndex,:))+StartPose(AgentIndex),numOfSite+1), 'r', 'LineWidth', 2);
ylabel("Space (5 m/site)");
yyaxis left
h = plot(timeWindowSet(1):seriesLenth+timeWindowSet(1)-1, 10*observabilityMetricSequence(:,1), 'k', ... 
    timeWindowSet(1):seriesLenth+timeWindowSet(1)-1, entropyOfSigmaGlobalSequence, 'b', ...
    timeWindowSet(1):seriesLenth+timeWindowSet(1)-1, conditionalEntropyOfSigmaGlobalGivenYSequence, 'g');
set(h, {'LineWidth'}, {2;2;2});
xlabel("Time (s/timestep)");
ylabel("Observability Metric");
legend('Metric \times 10', 'H(\Sigma)', 'H(\Sigma|Y)');
xlim([0 62]);
grid on







