% Rossler coupling for two cases:
% SSS(All oscillators have interaction through similar variables) and
% DDD(All oscillators have interaction through dissimilar varibales)


%



% SSS coupling
%Output comparison between full order and reduced order

figure;
clc;
clear;
t = 0:0.02:11;  % Time 
x30 = zeros(9,1);         %Initial conditions
C = zeros(9,1);               % C matrix
C(1) = 1;
[t,x3] = ode45( 'SSS_Coupling',t,x30); 
y3 = C'*x3';  


M = size(t,1);
y4 = zeros(M,3);

for i = 1:3
%set = [7,5,3];           % Order set is even
set = [8,6,4];          % Order set is odd


Ord = set(i);
%Initial conditions

[Ar Cr V] = Output_SSS(Ord);        % use function Output_matrix to get Cr and V 

x60 = V'*x30;
[t,x4] = ode45( 'SSS_Coupling_Reduced',t,x60); 
y4(:,i) = Cr*x4';                   % Output of reduced order systems
end




%%  Plots for the 3d trajectories
% plot3(x3(:,1),x3(:,4),x3(:,7));
% xlabel('x1','fontsize',24);
% ylabel('x2','fontsize',24);
% zlabel('x3','fontsize',24);
% grid on;
% 


plot(t,y3,t,y4(:,1),t,y4(:,2),t,y4(:,3));
%legend('9th order','8th order','6th order','4th order');
legend('9th order','7th order','5th order','3rd order');
xlabel('Simulation time(s)','fontsize',24);
ylabel('Response','fontsize',24);

grid on;





x1 =0.0001*ones(9,1);
mm(:,1) = x1;

%%
%  Start from 500 seconds later


tspan  = [0:0.5:550];   % time scale

M = size(tspan,2);
N = 9;
X = sym('x',[N 1]);



% x30 = 0.01*ones(9,1);     
C = zeros(9,1);
C(1) = 1;

%y = C'*x3';  

for i = 1:M

y = SSS_Coupling(tspan,x1);

mm(:,i+1) = 0.5*y+mm(:,i);

end




h0 = SSS_Coupling_Symbolic(tspan,X);
h = Lie_matrix(h0,X,N);           % Get the symbolic Lie derivative
% [t,y] = ode23( 'SSS_Coupling',tspan,x0);

OA = zeros(M,1);           % Full order observability

for j = M-100:M
    old = digits(5);
    Observ = vpa(subs(h,X,mm(:,j)));        % reduce the accuracy for calculation
    Eign = eig(Observ'*Observ);
    OA(j) = abs(min(Eign)/max(Eign));
end




OB = zeros(M,2);                  % Reduced order observability
%% -----------------------------------------------------------------------------------------
for i = 1:2
    
set = [5,3];            %Order set
N = 9;   
Ord = set(i);

X = sym('x',[Ord 1]);

% Initial conditions
[Ar Cr V] = Output_SSS(Ord);        % use function Output_matrix to get Cr 
x60 = V'*x1;

nn = zeros(Ord,2001);
nn(:,1) = x60;

for k = 1:M
z = SSS_Coupling_Reduced(tspan,x60);
nn(:,k+1) = 0.5*z+nn(:,k);
end

%[t,z] = ode45( 'SSS_Coupling_Reduced',tspan,x60);  
h0 = SSS_Coupling_Reduced_Symbolic(tspan,X);
h = Lie_matrix(h0,X,Ord);

%================================================================================

% Based on the observability matrix, we calculate the matrix eigenvalue
% ratio of min/max with time variation


for j = 1:M
    old = digits(5);
    Observ = vpa(subs(h,X,nn(:,j)));        % reduce the accuracy for calculation
    Eign = eig(Observ'*Observ);
    OB(j,i) = abs(min(Eign)/max(Eign));
end


end

semilogy(tspan,OA(:),tspan,OB(:,1),tspan,OB(:,2));
legend('Full order','5th order','3rd order');
xlabel('Simulation time','fontsize',24);
ylabel('Observability metric','fontsize',24);
grid on;
