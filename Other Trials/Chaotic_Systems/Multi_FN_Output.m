% % % Here I use two neuron system to do the simulation and compare the system
% % % observability after removing the states corresponding to fast dynamics
% % 
% % ========================================================================================================
% Response
Order_Series = [7,6,5,4];
t = 0:0.05:10;  % Time 
y = zeros(201,7);
FO = 50;


x0 = zeros(8,1);         %Initial conditions
C = zeros(1,8);               % C matrix
C(2) = 1;
[t,x_out1] = ode45('Multi_FN_Model',t,x0); 
y1 = C*x_out1';  

for i = 1:length(Order_Series)
    
    N = Order_Series(i);  
    x1 = zeros(8,1);         %Initial conditions
    C = zeros(1,8);               % C matrix
    
    C(2) = 1;
    if N == 4
        C(1) = 1;
        C(2) = 0;
    end
    x1(7) = N;
    [t,x_out2] = ode45('Reduced_Multi_FN',t,x1); 
    y = x_out2*C';  
    plot(t,y);
    hold on

end
hold on;
plot(t,y1);
legend('Full order','7th order','6th order', '5th order', '4th order');
xlabel('Simulation time','fontsize',24);
ylabel('System Output','fontsize',24);
grid on;

% 
% for i = 1:7
%     r(i) = 1-mean(abs(y1(10:end)'-y((10:end),i))./y1(10:end)');
% end
% plot(Order_Series,r);
% xlabel('System Order','fontsize',24);
% ylabel('Relative Accuracy','fontsize',24);
% grid on;
% 

%% Full order system ===========================================================================================================================================
% 
% 
% tspan  = [0:0.05:10];   % time scale
% % M = size(tspan,2);   %time length
% 
% M = size(tspan,2);
% FO = 4;
% N = FO;
% mm = zeros(3*N,M);    %Output for the full order system
% X = sym('x',[3*N 1]);    
% 
% x0 = zeros(FO*3,1);         %Initial conditions
% C = zeros(1,3*FO);               % C matrix
% mm(:,1) = x0;
% 
% % Get the state values at differnet time
% for s = 1:M
% y = Multi_Neuron(tspan(s),mm(:,s));
% mm(:,s+1) = 0.02*y+mm(:,s);
% end
% 
% h0 = Multi_Neuron_Sym(X);
% h = Neuron_Lie_matrix(h0,X,N);           % Get the symbolic Lie derivative
% 
% disp('step 1');
% Observ_Full = zeros(M,1);           % Full order observability
% % Substitute values into the symbols
% for j = 1:M
%      old = digits(5);
%      Observ = vpa(subs(h,X,mm(:,j)));        % reduce the accuracy for calculation
% %      f(X(:)) = h;
% %      FV = matlabFunction(h);
% %      Observ = FV(mm(:,j),mm(:,j),mm(:,j),mm(:,j),mm(:,j),mm(:,j));
%     Eign = eig(Observ'*Observ);
%     Observ_Full(j) = abs(min(Eign)/max(Eign));
% end
% 
% 
% %% Reduced order system
% % Order_Series = [18, 14, 12, 8, 6, 4, 3];
% Order_Series = [3,2,1];
% OB = zeros(3,401);
% disp('Full over'); 
% 
% for i = 1:length(Order_Series)
% disp('start');    
% FO = 20;
% N  = Order_Series(i);
% x1 = zeros(N*3,1);         %Initial conditions
% C = zeros(1,3*N);               % C matrix
% % C(3) = 1;
% % [t,x_out2] = ode45('Reduced_Multi_Neuron',t,x1); 
% % y(:,i) = x_out2*C';  
% 
% OR = FO-N;
% tspan  = [0:0.05:10];  % time scale
% M = size(tspan,2);         %time length
% X = sym('x',[3*N 1]);
% aa = zeros(3*N,M);    %Output for the full order system
% aa(:,1) = x1;
% 
% %  Get the state values at differnet time
% for ff = 1:M
%     y = Reduced_Multi_Neuron(tspan(i),aa(:,ff));
%     aa(:,ff+1) = 0.02*y+aa(:,ff);
% end
% 
% 
% h0 = Reduced_Multi_Neuron_Sym(X);
% h = Neuron_Lie_matrix(h0,X,N);           % Get the symbolic Lie derivative
% 
% Observ_Red = zeros(M,1);           % Reduced order observability
% 
% % Substitute values into the symbols
% for j = 1:M
%     old = digits(4);
%     Observ = vpa(subs(h,X,aa(:,j)));        % reduce the accuracy for calculation
%     Eign = eig(Observ'*Observ);
%     Observ_Red(j) = abs(min(Eign)/max(Eign));
% end
% OB(i,:) = Observ_Red';
% disp('over');  
% end 
% 
% semilogy(tspan,Observ_Full,tspan,Observ_Red);
% legend('Full order','Reduced Order');
% xlabel('Simulation time','fontsize',24);
% ylabel('Observability metric','fontsize',24);
% grid on;
% 
% 
% % 
% % [hAx,hLine1,hLine2] = plotyy(t,dd,t,b);
% % 
% %  xlabel('System orders','fontsize',24)
% % 
% %     ylabel(hAx(1),'Normalized accurate measurement','fontsize',24) % left y-axis
% %     ylabel(hAx(2),'Normalized Relative Observability','fontsize',24) % right y-axis
% %     grid on;  
% % 
% % 
% 
% 
% 
% 
% 
% 
% 
% 
