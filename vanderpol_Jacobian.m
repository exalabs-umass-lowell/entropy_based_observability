% Observability for x1/x2 by using the Jacobian matrix

% Q1 = [0   1;
%       -2*x1*x2-1   1-x1^2];       

% Q2 = [-(2*x1*x2+1)  1-x1^2;
%       (2*x1*x2+1)*(x1^2-1)  (x1^2-1)^2-1-2*x1*x2]

% P = Q1'*Q1
% Q = Q2'*Q2
% Then based on the min and max of eigenvalues from P/Q to
% estimate the observability index



%% Initial conditions and parameters
x0 = -10;
x1 = -8;

% Establish the van der pol model
t  = 0:0.01:6000;   % time scale
[t,x] = ode45( 'vd', t, [x0 x1],mu); 



%% Calculate the observability index

N = length(t);
D = zeros(N,1);
for i = 1:N
    
    % observability index for x1
    Q1 = [0 1;-2*x(i,1)*x(i,2)-1 1-x(i,1)^2];
    P = Q1'*Q1;
    A(i) = max(eig(P));
    B(i) = min(eig(P));
    C(i)=abs(B(i)/A(i));
    
    % observability index for x2
    Q2 = [-(2*x(i,1)*x(i,2)+1) 1-x(i,1)^2;
    (2*x(i,1)*x(i,2)+1)*(x(i,1)^2-1) (x(i,1)^2-1)^2-1-2*x(i,1)*x(i,2)];
    Q = Q2'*Q2;
    D(i) = max(eig(Q));
    E(i) = min(eig(Q));
    F(i)=abs(E(i)/D(i));

end



%% Plot the index along the trajectories

plot3(x(:,1),x(:,2),C);
xlabel('x1');
ylabel('x2');
zlabel('observability index');
title('variation of the x_1 observability with Jacobian matrix');
grid on;
figure;
plot3(x(:,1),x(:,2),F);
xlabel('x1');
ylabel('x2');
zlabel('observability index');
title('variation of the x_2 observability with Jacobian matrix');
grid on;