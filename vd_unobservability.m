% unobservability index for x2 by using the Lie Derivatives
% Os1 is invariant -> [1 0;0 1], thus the observability index of x1 is a constant. 
% Os2 = [0   1;
%       -2*mu*x1*x2  mu*(1-x1^2)]
% Q = Os2'*Os2 
% Then based on the reciprocal value of min eigenvalue to 
% estimate the unobservability index





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
    D(i)=abs(A(i)/B(i));
end

% plot the variations of observability index vs trajectories
figure;
plot3(x(:,1),x(:,2),D);
xlabel('x1', 'FontSize', 20);
ylabel('x2', 'FontSize', 20);
zlabel('observability index', 'FontSize', 20);
title('variation of the x_2 observability with Lie derivative', 'FontSize', 24);
grid on;









