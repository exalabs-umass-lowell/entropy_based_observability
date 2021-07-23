% Observability index for x2 by using the Lie Derivatives
% Os1 is invariant -> [1 0;0 1], thus the observability index of x1 is a constant. 
% Os2 = [0   1;
%       -2*mu*x1*x2  mu*(1-x1^2)]
% Q = Os2'*Os2 
% Then based on the min and max of eigenvalues from Q to
% estimate the observability index



%% Initial conditions and parameters
x0 = -10;
x1 = -8;


% Establish the van der pol model
t  = 0:0.01:6000;   % time scale
xa = [x0 x1];
[t,x] = ode45('vd', t, xa, mu); 



%% Calculate the observability index
Os1 = [1 0;0 1];
N = length(t);
A = zeros(N,1);
B = zeros(N,1);
D = zeros(N,1);
for i = 1:N
    Os2 = [0 1;-2*mu*x(i,1)*x(i,2)-1 mu*(1-x(i,1)^2)];  
    Q = Os2'*Os2;
    A(i) = max(eig(Q));
    B(i) = min(eig(Q));
    D(i)=abs(B(i)/A(i));
end

% plot the variations of observability index vs trajectories

plot3(x(:,1),x(:,2),D)
xlabel('x1');
ylabel('x2');
zlabel('observability index');
title('variation of the x_2 observability with Lie derivative');
grid on;









