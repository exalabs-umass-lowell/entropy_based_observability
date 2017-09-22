% % % % % Here I use Hindmarsh-Rose neuron system to do the simulation and compare the system
% % % % % observability after removing the states corresponding to fast dynamics
% % % %
% % % % ========================================================================================================
% % Response
Order_Series =  [7,5,3];
t = 0:0.02:50;  % Time
Time_Length = size(t,2);
y = zeros(Time_Length,3);
FO = 50;


x0 = zeros(3*3,1);         %Initial conditions
C1 = zeros(1,3*3);               % C matrix
C1(3) = 1;
[t,x_out1] = ode45('Neuron_Three_HR',t,x0);
y1 = C1*x_out1';


C2 = zeros(1,3*3);               % C matrix
C2(6) = 1;
y2 = C2*x_out1';

C3 = zeros(1,3*3);               % C matrix
C3(9) = 1;
y3 = C3*x_out1';


for i = 1:length(Order_Series)
    disp('start:')
    i
    
    N = Order_Series(i);
    M = (9-N)/2+N;
    x1 = zeros(M,1);         %Initial conditions
    C = zeros(1,M);               % C matrix
    C(3) = 1;
    x_out2 = zeros(Time_Length,M);
    for j= 1:Time_Length-1
        x_out2(j+1,:) = x_out2(j)+0.0001*Reduced_Neuron_Three_HR(t(j),x_out2(j,:)');
    end
    
%   [t,x_out2] = ode23s('Reduced_Neuron_Three_HR',t,x1);
    y(:,i) = x_out2*C';
    disp('End')
    i
end

% plot(t,y(:,1),t,y(:,2));
%plot(t,y1,t,y(:,1),t,y(:,2),t,y(:,3));
%legend('Full order','7th order','5th order', '3rd order');
plot3(y1,y2,y3);
xlabel('State x');
ylabel('State y');
zlabel('State z');
grid on;


for i = 1:3
    r(i) = 1-mean(abs(y1(10:end)'-y((10:end),i))./y1(10:end)');
end
plot(Order_Series,r);
xlabel('System Order','fontsize',24);
ylabel('Relative Accuracy','fontsize',24);
grid on;


%% Full order system ===========================================================================================================================================


tspan  = [0:0.05:20];   % time scale
% M = size(tspan,2);   %time length

M = size(tspan,2);
FO = 3;
N = FO;
mm = zeros(3*N,M);    %Output for the full order system
X = sym('x',[3*N 1]);

x0 = zeros(FO*3,1);         %Initial conditions
C = zeros(1,3*FO);               % C matrix
mm(:,1) = x0;

spmd
% Get the state values at differnet time
for s = 1:M
    y = Neuron_Three_HR(tspan(s),mm(:,s));
    mm(:,s+1) = 0.05*y+mm(:,s);
end

h0 = Neuron_Three_HR(0,X);
h = Neuron_Lie_matrix(h0,X,3*N);           % Get the symbolic Lie derivative

disp('Full Start');
Observ_Full = zeros(M,1);           % Full order observability
% Substitute values into the symbols
for j = 1:M
    old = digits(5);
    Observ = vpa(subs(h,X,mm(:,j)));        % reduce the accuracy for calculation
    %      f(X(:)) = h;
    %      FV = matlabFunction(h);
    %      Observ = FV(mm(:,j),mm(:,j),mm(:,j),mm(:,j),mm(:,j),mm(:,j));
    Eign = eig(Observ'*Observ);
    Observ_Full(j) = abs(min(Eign)/max(Eign));

end
end
disp('Full Over');

%% Reduced order system
% Order_Series = [18, 14, 12, 8, 6, 4, 3];
Order_Series = [7,5,3];
OB = zeros(3,401);

for i = 1:length(Order_Series)
    disp('start');
    i
    FO = 20;
    N  = Order_Series(i);
    x1 = zeros((9-N)/2+N,1);         %Initial conditions
    C = zeros(1,(9-N)/2+N);               % C matrix
    % C(3) = 1;
    % [t,x_out2] = ode45('Reduced_Multi_Neuron',t,x1);
    % y(:,i) = x_out2*C';
    
    OR = FO-N;
    tspan  = [0:0.05:50];  % time scale
    M = size(tspan,2);         %time length
    X = sym('x',[(9-N)/2+N 1]);
    aa = zeros((9-N)/2+N,M);    %Output for the full order system
    aa(:,1) = x1;
    
    %  Get the state values at differnet time
    for ff = 1:M
        y = Reduced_Neuron_Three_HR(tspan(i),aa(:,ff));
        aa(:,ff+1) = 0.02*y+aa(:,ff);
    end
    
    
    h0 = Reduced_Neuron_Three_HR(0,X);
    h = Neuron_Lie_matrix(h0,X,(9-N)/2+N);           % Get the symbolic Lie derivative
    
    Observ_Red = zeros(M,1);           % Reduced order observability
    
    % Substitute values into the symbols
    for j = 1:M
        old = digits(4);
        Observ = vpa(subs(h,X,aa(:,j)));        % reduce the accuracy for calculation
        Eign = eig(Observ'*Observ);
        Observ_Red(j) = abs(min(Eign)/max(Eign));
    end
    OB(i,:) = Observ_Red';
    disp('end')
    i
end

semilogy(tspan,Observ_Full,tspan,OB(1,:),tspan,OB(2,:),tspan,OB(3,:));
legend('9th Order','7th Order','5th Order','3rd Order');
xlabel('Simulation time','fontsize',24);
ylabel('Observability metric','fontsize',24);
grid on;


t = [9 7 5 3];
% r = mean(abs([y1(10:end)' y(10:end,1)' y(10:end,2)' y(10:end,3)']),1);
 % r = mean(abs([r1(2:80)' r2(2:80)' r3(2:80)' r4(2:80)' r5(2:80)']),1);
    dd = [1,r(1:3)'];
    dd = [1,(dd(2:4)-dd(4))/(1-dd(4))];
%     c = mean(OA,1);
%     b = ones(1,6)-(log(c/c(5))/log(c(1)/c(5)));
   
   [hAx,hLine1,hLine2] = plotyy(t,dd,t,[b(1:3),1]);

    xlabel('System orders','fontsize',24)

    ylabel(hAx(1),'Normalized accurate measurement','fontsize',24) % left y-axis
    ylabel(hAx(2),'Normalized Relative Observability','fontsize',24) % right y-axis
    grid on;  
%
% [hAx,hLine1,hLine2] = plotyy(t,dd,t,b);
%
%  xlabel('System orders','fontsize',24)
%
%     ylabel(hAx(1),'Normalized accurate measurement','fontsize',24) % left y-axis
%     ylabel(hAx(2),'Normalized Relative Observability','fontsize',24) % right y-axis
%     grid on;
%
%







%
