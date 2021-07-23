% Used to find the maximum values
%% We assume each small cell is square

t_start = 0;
t_end = 100;
num = 5000;
t = linspace(t_start,t_end,num);

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


v(1) = 0.4;
w(1) = -1.8;

   
%% van der pol limit cycle equations
%   .
%   v = w;
%   .
%   w   =  -0.08*(v^2-1)*w-v;


for(i = 1:1:num)
    v_dot(i) = w(i);
    w_dot(i) = -(v(i)^2-1)*w(i)-v(i);
    v(i+1) = v(i) + del_t*v_dot(i);
    w(i+1) = w(i) + del_t*w_dot(i);    
end
% figure(1)
%  plot(v,w,'r','Linewidth',2)
% hold on
% grid on;
% xlabel('x');
% ylabel('y');





%% Create an array of potential spatial scales to evaluate
Scales = [0.01 0.02 0.04 0.09 0.15 0.2 0.24 0.3];
% scale = linspace(0.1, 5, 5);
% max_time = 1000;


summ = zeros(length(Scales),1);
%%
for kk = 1:length(Scales)
    
   
% find the numbers in each axis
wid = ceil(7/Scales(kk));
len = ceil(7/Scales(kk));
M = zeros(wid,len);
xq = v;
yq = w;

    for i = 1:wid
        for j = 1:len
            xv = [-3.1+i*Scales(kk) -3+i*Scales(kk) -3+i*Scales(kk) -3.1+i*Scales(kk) -3.1+i*Scales(kk)];
            yv = [-3.1+0.5*j*Scales(kk) -3.1+0.5*j*Scales(kk) -3.0+0.5*j*Scales(kk) -3.0+0.5*j*Scales(kk) -3.1+0.5*j*Scales(kk)];
            [in, on] = inpolygon(xq,yq,xv,yv);
            M(i,j) = numel(xq(in));
        end
    end


% All the cells with data points inside are set as events    
[rows, colums, values] = find(M);
P = values./sum(values);    %Probability of even happening

% Probability of each event
Pd = length(values);

% Try to sum the entropy over all the possible events
for i = 1:Pd
    summ(kk) = summ(kk)+P(i)*log(1/P(i));
end

% Calculate the Entropy over per observation
summ(kk) = summ(kk)/Pd;

end

% plot the Entropy under different scales
plot(Scales,summ);
xlabel('Scales in Observation');
ylabel('Entropy per Observation');
grid on;
box on;




%% We assume each small cell is a rectangle
% 
% Scales = [0.2 0.5 1 2 4 8 12];
% 
% 
% for kk = 1:length(Scales)
%     
%     
%     
% Total_Number = 16*60;
% M = zeros(16,60);
% xq = quants;
% yq = quants2;
% 
%     for i = 1:16
%         for j = 1:60
%             xv = [-0.9+0.1*i -0.8+0.1*i -0.8+0.1*i -0.9+0.1*i -0.9+0.1*i];
%             yv = [-3.1+0.1*j -3.1+0.1*j -3.0+0.1*j -3.0+0.1*j -3.1+0.1*j];
%             [in, on] = inpolygon(xq,yq,xv,yv);
%             M(i,j) = numel(xq(in));
%         end
%     end
% 
%     
% [rows, colums, values] = find(M);
% P = values./650;
% summ(kk) = 0;
% for i = 1:55
%     summ(kk) = summ(kk)+P(i)*log(1/P(i));
% end
% 
% end
% 
% 
% plot(Scales,summ);
