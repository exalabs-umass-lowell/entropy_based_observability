% Code for the order reduction of Rossler systems
% by using Quadratic method to reduced the system order from 3 to 2
% Comparing states trajectories under different initial conditions to check if there
% exists some phenomenons related to chaos behavior


% Initial conditions varies between x = [-20 20], y = [-4 4]




    tspan  = 0:0.2:50;   % time scale
    x0 = [1 1 1];
    C = [1 0 0];
    [t,x1] = ode15s( 'Rossler',tspan,x0); 
    y1 = C*x1';
    [V Cr] = Para_Rossler();     % Output the parameters V and Cr
    x10 = x0*V;                 %Initial conditions
    [t,x2] = ode15s( 'Quadratic_Rossler',tspan,x10'); 
    y2 = Cr*x2';
    plot(t,y1,'-',t,y2,'--');
    legend('Full order','3rd order');
    xlabel('t','fontsize',24);
    ylabel('X','fontsize',24);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 20);
    box on;


    
