%% Observability of Rossler system for x/y/z by using the Lie derivative matrix

% A = Qx'*Qx
% C = Qy'*Qy
% E = Qz'*Qz
% Then based on the min and max of eigenvalues from P/Q to
% estimate the observability index


function Rossler_Lie_metric;

clc;

% Number of variable and initial conditions:

nbvar=3; 
xini=ones(1,nbvar)/10;

% Time parameters:

trans=100;
tend=5000;
tstep=0.01;
x0 = [1 1 1];



integration(x0,trans,tend,tstep);


%====================================================================
%====================================================================

function output=integration(x0,trans,tend,tstep);

[t,x] = run(x0,trans,tend,tstep);

a = 0.2; 
b = 0.2; 
c = 5.7; 

set(figure(2),'Position', [400 400 500 300]);  
clf;

plot3(x(:,1),x(:,2),x(:,3));
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
box on;
figure;

%% Calculate the observability index

N = length(t);
AB = zeros(N,1);
CD = zeros(N,1);
EF = zeros(N,1);

%% calculation of phase portraits /numerical solution/:
for i=1:N

        Qx = [ 1 0 0;0 -1 -1;-1-x(i,3) -a c-x(i,1)];
        Qy = [0 1 0;1 a 0;a a^2-1 -1];
        Qz = [0 0 1;x(i,3) 0 x(i,1)-c;b+2*x(i,3)*(x(i,1)-c) -x(i,3) (x(i,1)-c)^2-x(i,2)-2*x(i,3)];
        
        A(i) = min(eig(Qx'*Qx));
        B(i) = max(eig(Qx'*Qx));
        AB(i) = A(i)/B(i);
        C(i) = min(eig(Qy'*Qy));
        D(i) = max(eig(Qy'*Qy));
        CD(i) = C(i)/D(i);
        E(i) = min(eig(Qz'*Qz));
        F(i) = max(eig(Qz'*Qz));
        EF(i) = E(i)/F(i);
        
end



%% Plot the index along the trajectories

plot3(x(:,1),x(:,2),AB);
xlabel('x');
ylabel('y');
zlabel('observability index');
title('variation of the x observability with Lie derivative');
grid on;
figure;
plot3(x(:,1),x(:,2),CD);
xlabel('x');
ylabel('y');
zlabel('observability index');
title('variation of the y observability with Lie derivative');
grid on;
figure;
plot3(x(:,1),x(:,2),EF);
xlabel('x');
ylabel('y');
zlabel('observability index');
title('variation of the z observability with Lie derivative');
grid on;





% ===================================================================
% ===================================================================

function [t,x]=run(x0,trans,tend,tstep)

ttrans = [0:tstep:trans];
tspan = [0:tstep:tend];

option = [];

if trans > 0 
    [t x] = ode45(@dxdt,ttrans,x0,option);
    x0=x(end,:);
end

[t x] = ode45(@dxdt,tspan,x0,option);



% ===================================================================
% ===================================================================

function y = dxdt(t,x)

% parameters

a=0.2;
b=0.2;
c=5.7;

% equations

y = [
     -x(2)-x(3)+1;            % dx/dt
     x(1)+a*x(2);           % dy/dt
     b+x(3)*(x(1)-c);       % dz/dt
] ; 

