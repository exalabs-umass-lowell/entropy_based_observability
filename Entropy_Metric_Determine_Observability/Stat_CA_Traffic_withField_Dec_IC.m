% This is a function used to update the states of 1D CA traffic model

% Inputs cover:  intial configuration states and interaction coefficients J
% external coefficient B
% numSite ---  number of sites in the 1D dimension

function [Config] = Stat_CA_Traffic_withField_Dec_IC(J, B, scaling_Reps, numSite, timeLen, Config_init)

scale = 2;

timeLen = timeLen/scale^scaling_Reps;  % scale the time of simulation
GridNum = numSite/scale^scaling_Reps;  % scale the number of sites




% Each time scaled by 2, continuing scaling_Reps times of scaling
while(scaling_Reps >= 1)
    
    B = 1/2*log((exp(2*J+2*B)+exp(-2*J))/(exp(-2*J)+exp(2*J-2*B)));
    J = 1/4*log((exp(2*J+2*B)+exp(-2*J))*(exp(2*J-2*B)+exp(-2*J))/(exp(2*B)+2+exp(-2*B)));
    scaling_Reps = scaling_Reps-1;
end



% Scaling over the sensor range and road length

time_series = zeros(GridNum, timeLen);
Config = -ones(GridNum,timeLen);
Config(:,1) = Config_init;


for t = 1:timeLen-1
    i = 1;
    while(i<=GridNum)
        
        % In case get to the boundary
        ind = i;
        ind2 = i+1;
        if ind == GridNum
            ind2 = 1;
        end
        
        % Metroplis Algorithm to update the sites
        H = -J*Config(ind,t)*Config(ind2,t)-B*Config(ind,t);
        P_tr = exp(-H);        % probability based on Hamiltonian
        Rand_Val = rand;             % produce a random number
        
         % make sure the vehicle are moving forward/backward
         
        if Config(ind, t) == 1 && Config(ind,t)==-Config(ind2,t) && Rand_Val <= P_tr       
            Config(ind,t+1) = Config(ind2,t);
            Config(ind2,t+1) = Config(ind,t);
            i = i+1;
        else
            Config(ind,t+1) = Config(ind,t);
        end
        
        i = i+1;
    end
end




