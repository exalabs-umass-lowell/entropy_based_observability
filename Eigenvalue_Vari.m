%% In this m file, I try to observe the variations of eigenvalues for the A matrix by using linearization method








%% Full order system maximum eigenvalue variations

% a = 0.1;
% b = 0.1;
% c = 14;
% k = 0.176;
% 
% tspan  = 0:0.05:5;   % time scale
% M = size(tspan,2);
% N = 9;
% 
% % Calculate the Lie derivative matrix
% 
% 
% m = 0.005*ones(N,1);
% x30 = 0.01*ones(9,1);     
% C = zeros(9,1);
% C(1) = 1;
% [t,y] = ode45( 'SSS_Coupling',tspan,x30); 
% %y = C'*x3';  
%    
% % 
% % y = [-y1-z1+g1;           % Sinusoidal input on the first state
% %      x1+y1;
% %      b+z1*(x1-c);
% % 
% %     -y2-z2;
% %     x2+y2+g2;
% %     b+z2*(x2-c);
% % 
% % 
% %     -y3-z3;
% %     x3+y3;
% %     b+z3*(x3-c)+g3];
% % 
% 
% Eig = zeros(M,1);
% 
% for i = 1:M
%     x = y(i,:);
% A = [-k -1 -1 0 k 0 0 0 0;
%     1 1 0 0 0 0 0 0 0;
%     x(3) 0 x(1)-c 0 0 0 0 0 0;
%     0 0 0 0 -1 -1 0 0 0;
%     0 k 0 1 1-2*k 0 0 k 0;
%     0 0 0 x(6) 0 x(4)-c 0 0 0;
%     0 0 0 0 0 0 0 -1 -1;
%     0 0 0 0 0 0 1 1 0;
%     0 0 0 0 0 k x(9) 0 x(7)-k-c];
% 
% Eig(i) = max(real(eig(A)));
% 
% end
% 
% 
% plot(t,Eig);
% legend('Maxmum Eigenvalues');
% xlabel('Simulation time','fontsize',24);
% ylabel('Observability metric','fontsize',24);
% 




%%  Reduced order system maximum eigenvalue variations

a = 0.1;
b = 0.1;
c = 14;
k = 0.176;

tspan  = 0:0.05:5;   % time scale
M = size(tspan,2);
N = 9;

% Calculate the Lie derivative matrix


m = 0.005*ones(N,1);
x30 = 0.01*ones(9,1);     
%y = C'*x3';  


Eig = zeros(M,6);



for j = 1:6
    
    set = [7,6,5,4,3,2];            %Order set
    N = 9;   
    Ord = set(j);
    
    for i = 1:M

        X = sym('x',[Ord 1]);


        % Initial conditions
        [Ar Cr V] = Output_SSS(Ord);        % use function Output_matrix to get Cr 
        x60 = V'*x30;
        [t,z] = ode45( 'SSS_Coupling_Reduced',tspan,x60); 
        %z = Cr*x4';  

        x = z(i,:);

        h0 = SSS_Coupling_Reduced_Symbolic(tspan,X);
        h = jacobian(h0,X);

        %================================================================================

        M = size(t,1);

        % Based on the observability matrix, we calculate the matrix eigenvalue
        % ratio of min/max with time variation

        old = digits(5);
        Observ = vpa(subs(h,X,z(i,:)'));        % reduce the accuracy for calculation   
        Eig(i,j) = max(real(eig(Observ)));
       
    
    end

end




figure;
plot(t,Eig(:,1),t,Eig(:,2),t,Eig(:,3),t,Eig(:,4),t,Eig(:,5),t,Eig(:,6));
legend('7th Order', '6th Order', '5th Order', '4th Order', '3rd Order', '2nd Order');
xlabel('Simulation time','fontsize',24);
ylabel('Maximum Eigenvalues','fontsize',24);








