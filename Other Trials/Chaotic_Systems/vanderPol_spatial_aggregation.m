clear all
clc
close all

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

% Initial conditions
v(1) = 0.0;
w(1) = -0.1424;

%% van der pol limit cycle equations
%   .
%   x   =   x  - 0.3 x^3  - y  + u
%   .
%   y   =  0.08 ( x + 0.7 - 0.8 y)

for(i = 1:1:num)
    v_dot(i) = v(i) - (1/3)*(v(i))^3 - w(i) + I;
    w_dot(i) = 0.08*(v(i) + 0.7 - 0.8*w(i));
    v(i+1) = v(i) + del_t*v_dot(i);
    w(i+1) = w(i) + del_t*w_dot(i);    
end
figure(1)
plot(v,w,'r','Linewidth',2)
hold on;
grid on;

%% Spatial aggregation code
% Create an array of potential spatial scales to evaluate
scale = [0.2, 1, 5];
% scale = linspace(0.1, 5, 5);
max_time = 1000;

% for(x = 1:1:length(scale))
for(x = 1:1:length(scale))
%     scale(x) = 0.1;
    for(step = max_time+1:1:length(v)-max_time-1)
        v_arr = v(step-max_time:step+max_time);
        w_arr = w(step-max_time:step+max_time);
        
        r(step-max_time) = mean(v_arr(abs(v_arr - v(step)) < scale(x)));
        s(step-max_time) = mean(w_arr(abs(w_arr - w(step)) < scale(x)));
    end
    plot(r,s,'b.','Linewidth',2)
    pause(1)
end


plot(v,w,'r','Linewidth',2)