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
numOfAgent = 72;
interactionCoeff = 4.0;
externalFieldCoeff = 1.5;
lenOfTime = 64;

% Parameters for information quantification
influentialRange = 3;
timeWindowSet = [8, 12, 16];
AgentIndex = 3;


%% This function calculates the probability of the local states through Monte Carlo Simulations

CAT = data_process_util;
[Probability_Y_Given_SigmaM, interactionCoefficientVector, localStatesMatrix] = CAT.StatePredictionFromConditionalProbability(influentialRange, interactionCoeff, externalFieldCoeff);
% Calculate the time-space diagram
[Config, dataOfSpatialTemporal, StartPose, p_sequence] = CAT.StatisticalMechanicsBasedTraffic(numOfSite, numOfAgent, interactionCoefficientVector, externalFieldCoeff, lenOfTime, influentialRange);
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
        entropyOfSigmaGlobalAtOriginalScale = CAT.SigmaGlobalCalculation(localStatesMatrix, numOfSite, numOfAgent);
        [mutualInformation, observabilityMetric]  = CAT.observabilityQuantification(localStatesMatrix, Probability_SigmaM_Given_Y, numOfSite, numOfAgent, ProbabilityOfFreeFlowFromData, entropyOfSigmaGlobalAtOriginalScale);
        observabilityMetricSequence(t-t0+1, index) = observabilityMetric;
        mutualInformationSequence(t-t0+1, index)= mutualInformation;
    end
end

figure

subplot(1,2,1);
imagesc(Config);
set(gca,'YDir','normal');
set(gca,'FontSize',14)
colormap(gray);
xlabel("Time (s/site)");
ylabel("Space (5 m/site)");

hold on;

% test of motion data
AgentIndex = 3;
plot(1:lenOfTime, rem(cumsum(dataOfSpatialTemporal(AgentIndex,:))+StartPose(AgentIndex),numOfSite+1), 'r', 'LineWidth', 3);
xlabel("Time (s/timestep)");
ylabel("Space (5 m/site)");

% test of motion data
subplot(1,2,2);
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
xlim([0 lenOfTime]);
ylim([0 1]);
grid on


figure;
plot(1:lenOfTime, p_sequence, 'r.');
xlabel("Time (s/timestep)");
ylabel("Probability for exchange dynamics");
xlim([0 lenOfTime]);
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
xlim([0 lenOfTime]);
grid on



%%  Visualize on the variation of energy for the global systems

[n, m] = size(Config);
EnergyVar_Record = zeros(m,1);
for index = 1:m
    EnergyLevel = 0;
    for loc = 1:n-1
        EnergyLevel = EnergyLevel+Config(loc, index)*Config(loc+1, index);
    end
        EnergyVar_Record(index) = EnergyLevel;
end
plot(1:m, EnergyVar_Record, 'b', 'LineWidth', 2);

legend('Energy Variation', 'One observed trajectory');



