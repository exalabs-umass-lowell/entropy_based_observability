% We try to estimate the information from observation over the all the
% information on system states(condition entropy)

% Initial conditions
x0 = 0;
x1 = -2;


% Establish the van der pol model
t  = 0:0.05:1000;   % time scale
xa = [x0 x1];
[t,x] = ode45('vd', t, xa); 


% System trajectories for Vander Pol system
% 
% plot(x(:,1),x(:,2));
% xlabel('Simulation Time(s)');
% ylabel('State');
% 




%% Calculate the conditional entropy for state X
%=================================================================================================
N = length(x(:,1));







% Observation on X state
figure;
x_interval = -2:0.001:2;    % state interval
pd = fitdist(x(:,1),'kernel');
y = pdf(pd,x_interval);
plot(x_interval,y);
xlabel('State X');
ylabel('Probability');
grid on;






%% Calculate the conditional entropy for state X
%=================================================================================================



% Observation on X state
figure;
x_interval = -2:0.001:2;     % state interval
pd = fitdist(x(:,2),'kernel');
y = pdf(pd,x_interval);
plot(x_interval,y);
xlabel('State Y');
ylabel('Probability');
grid on;











X = x;
GMModel = fitgmdist(X,2);
figure
y = [zeros(1000,1);ones(1000,1)];
h = gscatter(X(:,1),X(:,2),y);
hold on
ezcontour(@(x1,x2)pdf(GMModel,[x1 x2]),get(gca,{'XLim','YLim'}))
title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
legend(h,'Model 0','Model1')
hold off



for j = 1:3
    subplot(2,2,j)
    h1 = gscatter(score(:,1),score(:,2),species);
    h = gca;
    hold on
    ezcontour(@(x1,x2)pdf(GMModels{j},[x1 x2]),...
        [h.XLim h.YLim],100)
    title(sprintf('GM Model - %i Component(s)',j));
    xlabel('1st principal component');
    ylabel('2nd principal component');
    if(j ~= 3)
        legend off;
    end
    hold off
end
g = legend(h1);
g.Position = [0.7 0.25 0.1 0.1];









