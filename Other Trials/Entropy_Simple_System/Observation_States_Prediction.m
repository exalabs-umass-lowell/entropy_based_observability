% We try to simulate the system dynamics of van der Pol and calculate the
% probability distribution

% Initial conditions
x0 = 0;
x1 = -2;


% Establish the van der pol model
t  = 0:0.05:1000;   % time scale
xa = [x0 x1];
[t,x] = ode45('vd', t, xa); 


%% System trajectories for Vander Pol system

plot(x(:,1),x(:,2));
xlabel('Simulation Time(s)');
ylabel('State');



%% Find the corresponding function to fit the distribution based on the observation from x/y direction over the system dynamics


% Observation on X state
figure;
x_interval = -2:0.001:2;    % state interval
pd = fitdist(x(:,1),'kernel');
y = pdf(pd,x_interval);
plot(x_interval,y);
xlabel('State X');
ylabel('Probability');
grid on;



% Observation on X state
figure;
x_interval = -2:0.001:2;     % state interval
pd = fitdist(x(:,2),'kernel');
y = pdf(pd,x_interval);
plot(x_interval,y);
xlabel('State Y');
ylabel('Probability');
grid on;




%%  Plot the bar distribution based on the observation from x/y direction over the system dynamics

figure;
y = x(:,1);
numIntervals = 100;     % number of observations on the system dynamics
intervalWidth = (max(y) - min(y))/numIntervals;
z = -1:intervalWidth:1;

%Next, use the HISTC function to find the frequency of each data range "z" in the given data set "y". This function returns the histogram count for a data set and range.
ncount = histc(y,z);

%Calculate the relative frequency of each data range by dividing the frequency by the total number of data points:
relativefreq = ncount/length(y);

%Finally plot the relative frequency versus the data ranges as a bar chart. On this chart, the bars will be adjoining, and the tick marks on the x-axis will label the extent of each bar's data range.
bar(z-intervalWidth/2, relativefreq,1);
xlabel('State X');
ylabel('Probability');
grid on;
xlim([min(z) max(z)]);



% State Y is similar to X above
figure;
y = x(:,2);

intervalWidth = (max(y) - min(y))/numIntervals;
z = -3:intervalWidth:3;

ncount = histc(y,z);
relativefreq = ncount/length(y);
bar(z-intervalWidth/2, relativefreq,1);
xlabel('State Y');
ylabel('Probability');
grid on;
xlim([min(z) max(z)]);
