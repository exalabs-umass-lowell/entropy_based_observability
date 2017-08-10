% We try to estimate the information from observation over the all the
% information on system states(condition entropy)

% Initial conditions
x0 = 0;
x1 = -2;


% Establish the van der pol model
t  = 0:0.05:200;   % time scale
xa = [x0 x1];
[t,x] = ode45('vd', t, xa); 


% System trajectories for Vander Pol system
% 
%plot(x(:,1),x(:,2));
% xlabel('Simulation Time(s)');
% ylabel('State');
%



%% Calculate the conditional entropy for state X
%=================================================================================================
N = length(x(:,1));
numIntervals_Group = [4 8 10 20 30 50 100 150 240 300 500];

y = x(:,1);

% number of observations on the system dynamics
for NI = 1:length(numIntervals_Group)

numIntervals_X = numIntervals_Group(NI);
    
interval_Length = (max(y) - min(y))/numIntervals_X;
State_X = -1:interval_Length:1;



%Next, use the HISTC function to find the frequency of each data range "z" in the given data set "y". This function returns the histogram count for a data set and range.
ncount = histc(y,State_X);
%ncount2 = histc(y,z);
%Calculate the relative frequency of each data range by dividing the frequency by the total number of data points:
relativefreq_X = ncount/length(y);


%% Calculate the conditional entropy for state Y
%=================================================================================================


% State Y is similar to X above

y = x(:,2);
numIntervals_Y = 3*numIntervals_X;
interval_Width = (max(y) - min(y))/numIntervals_Y;
State_Y = -3:interval_Width:3;

ncount = histc(y,State_Y);
relativefreq_Y = ncount/length(y);


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
        %[rows, colums, values] = find(Dot_Spread);

        if isempty(values)
            Prob_Value(i,j) = 0;
        else
            Prob_Value(i,j) = Dot_Spread(i,j)/L_X;    %Probability of event happening
        end
    end

end


% Normalize the probability distribution

Prob_Value = Prob_Value/sum(sum(Prob_Value));




%% Calculate the conditinal entropy along State X, X is known to predict the system states

% relativefreq_X = relativefreq_X/sum(relativefreq_X);
% relativefreq_X = relativefreq_X/sum(relativefreq_Y);

EY_X = 0;

for i = 1:L_X-1
    for j = 1:L_Y-1
        if relativefreq_X(i)~=0 && Prob_Value(i,j)~= 0
            EY_X = EY_X + Prob_Value(i,j)*log(Prob_Value(i,j)/relativefreq_X(i));
        end
    end 
    
end

% Variation in different scales when known x
EY(NI) = EY_X;
%% Calculate the conditinal entropy along State Y, Y is known to predict the system states

EX_Y = 0;

for i = 1:L_X-1
    for j = 1:L_Y-1
        if relativefreq_Y(j)~= 0 && Prob_Value(i,j)~= 0
            EX_Y = EX_Y + Prob_Value(i,j)*log(Prob_Value(i,j)/relativefreq_Y(j));
        end
        
    end 
    
end


% Variation in different scales when known y
EX(NI) = EX_Y;

end

plot(numIntervals_Group,EX,numIntervals_Group,EY);
xlabel('Scales');
ylabel('Condition Entropy');
legend('Known X', 'Known Y');
grid on;



