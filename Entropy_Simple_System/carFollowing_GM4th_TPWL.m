% Here we use the model of GM4th model to do the simulation on traffic inside a
% ring road, acceleration a_n(t) = c*v_n(t) * Delta_v(t-T)/Delta_x(t-T)

% Simulation of 20 vehicles using the GM4th model for the closed ring road



% Initial conditions of the traffic
Radius = 12;
RingTrack_length = 2*pi*Radius;
Vehicle_Num = 20;

safeTimeGap = 1.6;
AccExponent = 4;

v0 = 5;  %velocity
s0 = 2;  %position

% Initial states:   n positions and n velocities
x0 = zeros(2*Vehicle_Num,1);

%%   Initial conditions for the traffic
    

% from 1 to vehicle_num, correspoding to the position of each vehicle

for i = 1:Vehicle_Num/2
    x0(i) = 10+(s0+2)*(Vehicle_Num-i);
end


for i = (Vehicle_Num/2)+1:Vehicle_Num
    x0(i) = (s0+1)*(Vehicle_Num-i);
end

% 
% from vehicle_num+1 to 2*vehicle_num, correspoding to the velocity of each vehicle

for i = Vehicle_Num+1:Vehicle_Num+(Vehicle_Num/2)
    x0(i) = 4;
end

for i = Vehicle_Num+(Vehicle_Num/2)+1:2*Vehicle_Num
    x0(i) = 12;
end


%Acceleration_vector = zeros(2*Vehicle_Num,1);
s = zeros(Vehicle_Num,1);

incrTime = 0.02; % Time step for numerical integration
% yout = zeros(6000,2*Vehicle_Num);
% yout(1,:) = x0(1:2*Vehicle_Num)';
% Acceleration_vector = GM_4th(yout(1,:),Vehicle_Num,RingTrack_length);
% yout(2,:) = yout(1,:)+Acceleration_vector'*incrTime;
% Acceleration_vector = GM_4th(yout(2,:),Vehicle_Num,RingTrack_length);
% yout(3,:) = yout(2,:)+Acceleration_vector'*incrTime;



tspan = linspace(0,120,6000);     % Two minute time


%% Select N reference points based on the K-means algorithm
% Here the K-means algorithm is from "Statistical and Machine Learning
% Toolbox"


Cluster_NumSeries = linspace(2,20,19);  % the series used to find the optimal cluster size
Sum_SquareError = zeros(19,1);
for ClusNum = 1:length(Cluster_NumSeries)

    [idx, C, sumd] = kmeans(yout,ClusNum,'MaxIter',3000);
    Sum_SquareError(ClusNum) = sum(sumd);

end

plot(Cluster_NumSeries,Sum_SquareError);
xlabel('Number of clusters', 'FontSize', 20);
ylabel('Within Clusters Sum of Squares', 'FontSize', 20);
grid on;
box on;
%%==================================================================================

% Based on the above section, we determined the number of reference points
% to be 6

% figure;
ClusterNum = 6;
[idx, Ref_x, sumd] = kmeans(yout,ClusterNum,'MaxIter',3000);   % Ref_x is the states at each reference points/ centroid for the cluster

% f(x) at all the reference points form the vector
Ref_FX = zeros(ClusterNum,2*Vehicle_Num); % Initialization

for i = 1:ClusterNum
    Ref_FX(i,:) = GM_4th(Ref_x(i,:),Vehicle_Num, RingTrack_length);
end

%%
%==================================================================================================================


% Here we try to use the 6 reference points to approximate the acceleration
% values for each trajectory points

Time_length = size(tspan,2);

%  f(x) vector, used to approximate the real f(x)
FX_Approx = zeros(ClusterNum, 2*Vehicle_Num);
x = zeros(Time_length, 2*Vehicle_Num);   %States series of the traffic model
x(1,:) = x0;   %Initial states

for e = 1:Time_length-1
        
%Calculate the distances between points and the reference ones
beta = 25;   % a constant
Distance = zeros(ClusterNum,1);
for i = 1:ClusterNum
    Distance(i) = sqrt(sum((x(e,:) - Ref_x(i,:)).^2));    
end

% Based on the distance, we try to calculate the weights for each reference
% point
m = min(Distance);
Reference_Weights = exp(-beta*Distance(:)./m);  

% Normalization over the weights
WeightSum = sum(Reference_Weights);
Reference_Weights = Reference_Weights(:)/WeightSum; 


% dx/dt = f(x0) + A*(x-x0);           A = [A1; A2]
A1 =  [zeros(Vehicle_Num) eye(Vehicle_Num)];
A2 = zeros(Vehicle_Num,4);


for i = 1:ClusterNum
    A21 = PointsSelected_Linearization(Ref_x(i,:),RingTrack_length);
    FX_Approx(i,1:Vehicle_Num) = x(e,Vehicle_Num+1:2*Vehicle_Num);

    % A = [A1;A2];
    for s = 1+Vehicle_Num:Vehicle_Num*2-1
        
        %  Taylor series expansion, taking the constant and linear terms
        FX_Approx(i,s) = Ref_FX(i,s) + A21(s-Vehicle_Num,:)*[x(e,s-Vehicle_Num)-Ref_x(i,s-Vehicle_Num);...
        x(e,s-Vehicle_Num)-Ref_x(i,s-Vehicle_Num);x(e,s)-Ref_x(i,s);x(e,s+1)-Ref_x(i,s+1)];

    end

end

% State at next time point is  x(i) = x(i-1) + dx * dt, here dx is equal to
% combination of dx at reference points 

x(e+1,:) = x(e,:) + Reference_Weights'*FX_Approx(:,:).*incrTime;

end


hold all;

for i = 1:Vehicle_Num
    plot(tspan, x(:,i));
    hold on;
end

xlabel('Time', 'FontSize', 20);
ylabel('Position', 'FontSize', 20);
%title('Simulation of GM 4th in closed boundary road', 'FontSize', 24);

grid on;
box on;






