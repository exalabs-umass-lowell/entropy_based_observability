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
numOfAgent = 16;
interactionCoeff = 4.0;
externalFieldCoeff = 1.5;
lenOfTime = 64;
scalingConst = 2;



numOfSite_2 = numOfSite/scalingConst;
numOfSite_3 = numOfSite_2/scalingConst;
numOfAgent
numOfAgent




[Config, dataOfSpatialTemporal] = CAT.StatisticalMechanicsBasedTraffic(numOfSite, numOfAgent, interactionCoefficientVector, externalFieldCoeff, lenOfTime, influentialRange);



[Config_2, dataOfSpatialTemporal] = CAT.StatisticalMechanicsBasedTraffic(numOfSite, numOfAgent, interactionCoefficientVector, externalFieldCoeff, lenOfTime, influentialRange);


[Config_3, dataOfSpatialTemporal] = CAT.StatisticalMechanicsBasedTraffic(numOfSite, numOfAgent, interactionCoefficientVector, externalFieldCoeff, lenOfTime, influentialRange);




%%  Results visualization

subplot(3,1,1)
imagesc(Config);
set(gca,'YDir','normal');
colormap(gray);
xlabel("Time (s/site)");
ylabel("Space (5 m/site)"); 

subplot(3,1,2)
imagesc(Config_2);
set(gca,'YDir','normal');
colormap(gray);
xlabel("Time (s/site)");
ylabel("Space (5 m/site)"); 

subplot(3,1,3)
imagesc(Config_3);
set(gca,'YDir','normal');
colormap(gray);
xlabel("Time (s/site)");
ylabel("Space (5 m/site)"); 




%% Output the values/resutls

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
