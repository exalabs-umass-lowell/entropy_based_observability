% Code to get the output from sinusoidal input on x state

    
   tspan  = 0:0.01:50;   % time scale
   
   [t,x1] = ode15s( 'HyperRossler',tspan,[-10,-6,0,10]);
%    hold on;
%   [t,x2] = ode15s( 'HyperRossler_SingularPerturb',tspan,[4;105;1]);
  
   
%     plot(t,x1(:,1),'--',t,x2(:,1),'o');
%     legend('Full order','3rd order');
%     xlabel('t','fontsize',24);
%     ylabel('X','fontsize',24);
%     xt = get(gca, 'XTick');
%     set(gca, 'FontSize', 20);
%     grid on;
%     box on;

 plot3(x1(:,1),x1(:,2),x1(:,3));
 xlabel('x1','fontsize',24);
 ylabel('x2','fontsize',24);
 zlabel('x3','fontsize',24);
grid on;
