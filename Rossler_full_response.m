%% Observability of Rossler system for x/y/z by using the Lie derivative matrix

% A = Qx'*Qx
% C = Qy'*Qy
% E = Qz'*Qz
% Then based on the min and max of eigenvalues from P/Q to
% estimate the observability index


function Rossler_full_response;

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

a = 0.25; b = 0.2; c = 1; 

clf;
plot(t,x(:,1));
xlabel('t','fontsize',18);
ylabel('y','fontsize',18);
box on;
figure;




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

a = 0.25; b = 0.2; c = 1;

% equations

y = [
     -x(2)-x(3)+1;            % dx/dt
     x(1)+a*x(2);           % dy/dt
     b+x(3)*(x(1)-c);       % dz/dt
] ; 

