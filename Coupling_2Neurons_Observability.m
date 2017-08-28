% Here I use two neuron system to do the simulation and compare the system
% observability after removing the states corresponding to fast dynamics


t = 0:0.05:60;  % Time 
x0 = zeros(6,1);         %Initial conditions
C = zeros(1,6);               % C matrix
C(3) = 1;
[t,x_out1] = ode45('Neuron_System',t,x0); 
y1 = C*x_out1';  


x1 = zeros(4,1);         %Initial conditions
C = zeros(1,4);               % C matrix
C(2) = 1;
[t,x_out2] = ode45('ReducedNeuron_System',t,x1); 
y2 = C*x_out2';  

% plot(t,y1,t,y2);
% legend('Full order','Reduced order');
% xlabel('Simulation time','fontsize',24);
% ylabel('System Output','fontsize',24);
% grid on;
% 
% 


%% Full order system
tspan  = [0:0.05:5];   % time scale
M = size(tspan,2);   %time length
mm = zeros(6,M);    %Output for the full order system
M = size(tspan,2);
N = 6;
X = sym('x',[N 1]);    


%  Get the state values at differnet time
for i = 1:M
y = Neuron_System(tspan(i),mm(:,i));
mm(:,i+1) = 0.05*y+mm(:,i);
end


h0 = Neuron_System(tspan,X);
h = Lie_matrix(h0,X,N);           % Get the symbolic Lie derivative

Observ_Full = zeros(M,1);           % Full order observability
% Substitute values into the symbols
for j = 1:M
    old = digits(5);
    Observ = vpa(subs(h,X,mm(:,j)));        % reduce the accuracy for calculation
    Eign = eig(Observ'*Observ);
    Observ_Full(j) = abs(min(Eign)/max(Eign));
end


%% Reduced order system

tspan  = [0:0.05:5];   % time scale
M = size(tspan,2);         %time length
mm = zeros(4,M);             % Output for the reduced system
M = size(tspan,2);
N = 4;
X = sym('x',[N 1]);
 
%  Get the state values at differnet time
for i = 1:M
y = ReducedNeuron_System(tspan,mm(:,i));
mm(:,i+1) = 0.05*y+mm(:,i);
end


h0 = ReducedNeuron_System(tspan,X);
h = Lie_matrix_Reduced(h0,X,N);           % Get the symbolic Lie derivative

Observ_Red = zeros(M,1);           % Full order observability

% Substitute values into the symbols
for j = 1:M
    old = digits(5);
    Observ = vpa(subs(h,X,mm(:,j)));        % reduce the accuracy for calculation
    Eign = eig(Observ'*Observ);
    Observ_Red(j) = abs(min(Eign)/max(Eign));
end


semilogy(tspan,Observ_Full,tspan,Observ_Red);
legend('Full order','Reduced Order');
xlabel('Simulation time','fontsize',24);
ylabel('Observability metric','fontsize',24);
grid on;





