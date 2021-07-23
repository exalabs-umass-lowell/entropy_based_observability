%% Observability of Rossler system for x/y/z by using the Lie derivative matrix

% A = Qx'*Qx
% C = Qy'*Qy
% E = Qz'*Qz
% Then based on the min and max of eigenvalues from P/Q to
% estimate the observability index


function x = Rossler_orderReduction();

clc;

% Number of variable and initial conditions:

nbvar=3; 
xini=ones(1,nbvar)/10;

% Time parameters:

trans=100;
tend=3000;
tstep=0.1;
x0 = [1 1 1];



x = integration(x0,trans,tend,tstep);


%====================================================================
%====================================================================

function x=integration(x0,trans,tend,tstep);

[t,x] = run(x0,trans,tend,tstep);

a = 0.2; 
b = 0.2; 
c = 5.7; 



plot3(x(:,1),x(:,2),x(:,3));
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
box on;
grid on;





% ===================================================================
%% ===================================================================

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
     -x(2)-x(3);            % dx/dt
     x(1)+a*x(2);           % dy/dt
     b+x(3)*(x(1)-c);       % dz/dt
] ; 

