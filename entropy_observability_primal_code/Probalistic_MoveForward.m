% This is a function used to update the states of 1D CA traffic model


function Config = Probalistic_MoveForward(randomness, numSite, timeLen, Config_init)



GridNum = numSite;  

% Scaling over the sensor range and road length

time_series = zeros(GridNum, timeLen);
Config = ones(GridNum,timeLen);
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
        
        Rand_Val = rand;
        
         % make sure the vehicle are moving forward/backward
         
        if(randomness>0 && ind == 1)
            P_tr = 0.5;
        else 
            P_tr = 1.0;
        end
        
        if Config(ind,t) == 1 && Config(ind,t)==-Config(ind2,t) && Rand_Val <= P_tr       
            temp = Config(ind,t);
            Config(ind,t+1) = Config(ind2,t);
            Config(ind2,t+1) = temp;
            i = i+1;
            time_series(ind, t+1) = 1;
        else
            Config(ind,t+1) = Config(ind,t);
        end
        
        i = i+1;
    end
end




