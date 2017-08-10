% Observability of Rossler system for x/y/z by using the Lie Derivative

% A = Qx'*Qx
% C = Qy'*Qy
% E = Qz'*Qz
% Then based on the min and max of eigenvalues from A/C/E to
% estimate the observability index



%% Initial conditions and parameters
a = 0.2; 
b = 0.2; 
c = 5.7;    
Y0 = [1 1 0];
% Establish the van der pol model
t = 0:0.01:50;   % time scale
[t,x] = ode45( @rs, t, Y0); 



%% Calculate the observability index


N = length(t);
AB = zeros(N,1);
CD = zeros(N,1);
EF = zeros(N,1);

%% calculation of phase portraits /numerical solution/:
for i=2:N

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
xlabel('x1');
ylabel('x2');
zlabel('observability index');
title('variation of the x observability with Lie derivative');
grid on;
figure;
plot3(x(:,1),x(:,2),CD);
xlabel('x1');
ylabel('x2');
zlabel('Observability Index');
%title('variation of the y observability with Lie derivative');
grid on;
figure;
plot3(x(:,1),x(:,2),EF);
xlabel('x1');
ylabel('x2');
zlabel('Observability Index');
title('variation of the z observability with Lie derivative');
grid on;

figure;
plot3(x(:,1),x(:,2),x(:,3));
xlabel('x1');
ylabel('x2');
zlabel('x3');
grid on;