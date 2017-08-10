%% Initial conditions and parameters
a = 0.2; 
b = 0.2; 
c = 5.7;    
Y0 = [0 5 1];
% Establish the van der pol model
t = [0,120];   % time scale
[t,x] = ode45( @rs, t, Y0); 


%% From "?????" plane to observe the system

%=================================================================================================
N = length(x(:,1));
% numIntervals_Group = [4 8 10 20 30 50 100 150 240 300 500];
numIntervals_Group = 50;




% %  Mapping to X-Y plane 
% m = x(:,1);
% n = x(:,2);


%  Mapping to X-Y plane 
m = x(:,2);
n = x(:,3);
% 
% %  Mapping to X-Y plane 
% m = x(:,1);
% n = x(:,3);




% number of observations on the system dynamics
for NI = 1:length(numIntervals_Group)

    
numIntervals_X = numIntervals_Group(NI);
    
interval_Length = (max(m) - min(m))/numIntervals_X;
State_X = min(m):interval_Length:max(m);



%Next, use the HISTC function to find the frequency of each data range "z" in the given data set "y". This function returns the histogram count for a data set and range.
ncount = histc(m,State_X);
%ncount2 = histc(y,z);
%Calculate the relative frequency of each data range by dividing the frequency by the total number of data points:
relativefreq_X = ncount/length(m);


%% Calculate the conditional entropy for state Y
%=================================================================================================


% State Y is similar to X above


numIntervals_Y = (max(n) - min(n))/interval_Length;

State_Y = min(n):interval_Length:max(n);

ncount = histc(n,State_Y);
relativefreq_Y = ncount/length(n);


L_X = length(State_X);
L_Y = length(State_Y);
Dot_Spread = zeros(L_X,L_Y);
Prob_Value = zeros(L_X,L_Y);

xq = x(:,1);
yq = x(:,2);
      
T_Num = 0;



%% Calculate the probability distribution




k = 0;
for i = 1:L_X-1
    for j = 1:L_Y-1


        % Use inpolygon function to discretize them and count the number of spots
        % inside each cell

        xv = [State_X(i) State_X(i+1) State_X(i+1) State_X(i) State_X(i)];
        yv = [State_Y(j) State_Y(j) State_Y(j+1) State_Y(j+1) State_Y(j)];
        [in, on] = inpolygon(xq,yq,xv,yv);
        
        Dot_Spread(i,j) = numel(xq(in));
        %values = sum(in);
           
        % All the cells with data points inside are set as events
        [rows, colums, values] = find(Dot_Spread);

        if isempty(values)
            Prob_Value(i,j) = 0;
        else
            Prob_Value(i,j) = Dot_Spread(i,j)/L_X;    %Probability of event happening
        end
    end

end


% Normalize the probability distribution

Prob_Value = Prob_Value/sum(sum(Prob_Value));



end

figure;

mesh(State_X, State_Y, Prob_Value');
xlabel('State Y');
ylabel('State Z');
title('Y-Z plane');
grid on;
box on;


