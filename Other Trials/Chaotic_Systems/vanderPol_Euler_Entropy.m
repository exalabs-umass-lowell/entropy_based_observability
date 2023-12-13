
% In this code, I try to set the initial conditions of all the dots
% distributed uniformly in the square x = [0,1], y = [0,1];

clc;
clear all;
t_start = 0;
t_end = 1;
num = 5000;
t = linspace(t_start,t_end,num);


% Total number of points
N = 5000; 
%%
% van der Pol model parameters
a = 0.8;
b = 0.4;
c = 0.2;
I = 0.33;

% Array initialization
del_t = (t_end-t_start)/num;
v_dot = zeros(1,length(t));
v = v_dot;
w = v_dot;
w_dot = v_dot;


% Uniformly distributed initial conditions

x_0 = rand(1, N)-0.1;
y_0 = rand(1, N)-2.3;
minAllowableDistance = 0.05;
numberOfPoints = N;
% Initialize first point.
keeperX = x_0(1);
keeperY = y_0(1);
% Try dropping down more points.
counter = 2;
for k = 2 : numberOfPoints
	% Get a trial point.
	thisX = x_0(k);
	thisY = y_0(k);
	% See how far is is away from existing keeper points.
	distances = sqrt((thisX-keeperX).^2 + (thisY - keeperY).^2);
	minDistance = min(distances);
	if minDistance >= minAllowableDistance
		keeperX(counter) = thisX;
		keeperY(counter) = thisY;
		counter = counter + 1;
	end
end
plot(keeperX, keeperY, 'b*');

% grid on;
% hold on;
N = size(keeperX,2);
% 
%Scales = linspace(0.03, 4,10);
Scales = [0.3 0.6 1 1.5 2 2.5 2.8 3.2];
Scales = [0.5];
% v = zeros(100,num);
% w = zeros(100,num);

Entropy_Time = zeros(length(Scales),14);
    
%% van der pol limit cycle equations
%   .
%   v = w;
%   .
%   w   =  -0.08*(v^2-1)*w-v;
for(i = 1:1:num)
    
    
    for mm = 1:N
        
        v(mm,1) = keeperX(mm);
        w(mm,1) = keeperY(mm);
        
        v_dot(i) = w(mm,i);
        w_dot(i) = -(v(mm,i).^2-1)*w(mm,i)-v(mm,i);
        v(mm,i+1) = v(mm,i) + del_t*v_dot(i);
        w(mm,i+1) = w(mm,i) + del_t*w_dot(i);
        %scatter(v(:,i),w(:,i));
        %pause(0.1);
        
    end
    
end



% plot(v(5,:),w(5,:),'r','Linewidth',2);
% hold on;
% 
% for(i = 1:1:num)
%     
% plot(v(5,:),w(5,:),'r','Linewidth',2);
% hold on;
% scatter(v(:,i),w(:,i));
% pause(0.1);
% hold off
% end 

%%
for(i = 1:300:4201)
    
    summ = zeros(length(Scales),1);
    
    
    for kk = 1:length(Scales)
        
        
        % find the numbers in each axis
        wid = ceil(10/Scales(kk));
        len = ceil(10/Scales(kk));
        M = zeros(wid,len);
        xq = v(:,i);
        yq = w(:,i);
        
        
        % Use inpolygon function to discretize them and count the number of spots
        % inside each cell
        
        for ii = 1:wid
            for j = 1:len
                xv = [-5+ii*Scales(kk) -5+(ii+1)*Scales(kk) -5+(ii+1)*Scales(kk) -5+ii*Scales(kk) -5+ii*Scales(kk)];
                %yv = [-4+0.5*j*Scales(kk) -3.1+0.5*j*Scales(kk) -3.0+0.5*j*Scales(kk) -3.0+0.5*j*Scales(kk) -3.1+0.5*j*Scales(kk)];
                yv = [-5+j*Scales(kk) -5+j*Scales(kk) -5+(j+1)*Scales(kk) -5+(j+1)*Scales(kk) -5+j*Scales(kk)];
                [in, on] = inpolygon(xq,yq,xv,yv);
                M(ii,j) = numel(xq(in));
            end
        end
        
        
        % All the cells with data points inside are set as events
        [rows, colums, values] = find(M);
        P = values./sum(values);    %Probability of event happening
        
        %
        Pd = length(values);
        
        % Try to sum the entropy over all the possible events
        for s = 1:Pd
            summ(kk) = summ(kk)+P(s)*log(Scales(kk)^2*1/P(s));
        end

        
        % Calculate the Entropy over per observation
        summ(kk) = summ(kk)/Pd;
        
    end
    
    Entropy_Time(:,i) = summ;
    
end



Entropy = mean(Entropy_Time,2);

%plot the Entropy under different scales
plot(Scales,Entropy./Scales'.^2);
xlabel('Scales in Observation');
ylabel('Entropy per Observation');
grid on;
box on;

