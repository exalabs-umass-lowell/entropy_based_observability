%  Observability variation caused by order reduction for a single neuron system

%% Output of system comparison before and after order reduction


t = 0:0.02:30;  % Time
Time_Length = size(t,2);
y = zeros(Time_Length,3);

% x0 = [5,-4,0];        %Initial conditions
x0 = zeros(1,3);
C = zeros(1,3);               % C matrix
C(3) = 1;
[t,x_out1] = ode45('Single_Neuron',t,x0);
y1 = C*x_out1';


x0 = zeros(2,1);
x0(2) = 0;   %Initial conditions
C = [0,1];
[t,x_out2] = ode45('Reduced_Single_Neuron',t,x0);
y2 = C*x_out2';

plot(t,y1,t,y2);
legend('Full order','Reduced order');
xlabel('Simulation time','fontsize',24);
ylabel('System Output','fontsize',24);
grid on;
box on;


%% Observability comparison before and after order reduction
% 
% t = 0:0.02:30;  % Time
% Time_Length = size(t,2);
% y = zeros(Time_Length,3);
% 
% x0 = [1,-4,1];        %Initial conditions
% C = zeros(1,3);               % C matrix
% C(3) = 1;
% [t,x_out1] = ode45('Single_Neuron',t,x0);
% y1 = C*x_out1';
% 
% x0 = ones(2,1);         %Initial conditions
% C = [0,1];
% [t,x_out2] = ode45('Reduced_Single_Neuron',t,x0);
% y2 = C*x_out2';
% 
% 
% plot(t,y1,t,y2);
% legend('Full order','Reduced order');
% xlabel('Simulation time','fontsize',24);
% ylabel('System Output','fontsize',24);
% grid on;
% box on;
% 
% 
