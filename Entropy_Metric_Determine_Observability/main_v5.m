%% Build time-space diagrams

close all
clear all
clc

flag_plot = 0;
flag_plot_noCV = 0;
flag_plot_delay_profile = 1;
flag_plot_dela = 0;
flag_plot_N_t = 0;
figProp = 125;

count_param = 1;

%% Define parameters
qmax = 1800;        % Capacity (veh/hr)
qc = qmax;
vf = 90;            % free flow speed (km/hr)
kj = 110;           % jam density (veh/km)
params.kj = kj;
params.qc = qc;
params.vf = vf;

kc = qmax/vf;       % Critical density
params.kc = kc;
w = qc/(kj-kc);
params.w = w;
qj = 0;
params.qj = qj;

% Traffic state parameters
% ka = 0.825*kc;
ka = 0.9*kc;
qa = ka*vf;
params.ka = ka;

kg = 0.5*kc;
% kg = 0.7*kc;
ke = get_congested_density(kg*vf, params);
% Bottleneck capacity is reduced by 15% at state I/H
ki = 0.85*kg;
kh = get_congested_density(ki*vf, params);
kf = 0.25*kc; %TODO: Must be independent of bottleneck capacity
% kf = 0.2975*kc;
qf = kf*vf;

horizon_null = [];
horizon_event = [];

%% LOOP 1: Parameter being changed
for ka = 0.5*kc:(0.01*kc):kc
% for kf = 0.01*kc:0.01*(kc-ki):ki
% for kg = 1.5*kf:(0.05*(ka-kf)):ka
% for kg = 0.7*kc:(0.02*(ka-kf)):0.7*kc
    count_delay = 1;
    ka_data(count_param) = ka;
    qa = ka*vf;
    
%     kf_data(count_param) = kf;
%     
%     kg_data(count_param) = kg;
    ke = get_congested_density(kg*vf, params);
    % Bottleneck capacity is reduced by 15% at state I/H
    ki = 0.85*kg;
    kh = get_congested_density(ki*vf, params);
%     kf = 0.25*kc; %TODO: Must be independent of bottleneck capacity
    qf = kf*vf;

    
    ve = get_vel(ke, params);
    qe = ke*ve;
    qg = qe;
    vh = get_vel(kh, params);
    qh = kh*vh;
    qi = qh;
    
    vae = -(qa - ke*ve)/(ka-ke);
    vef = (qe-qf)/(ke-kf);
    vfh = (qh-qf)/(kh-kf);
    vah = -(qa-qh)/(ka-kh);   
    
    
    %% Define traffic states
    T(1).value = 0;
    T(1).label = 'T_0';
    T(1).descriptor = 'Initial time';
    
    x1(1).value = 5; %5;     % km
    x1(1).descriptor = 'Initial position of CV1';
    x2(1).value = 5.1; %5.1;   % x2 > x1 (always)
    x2(1).descriptor = 'Initial position of CV2';
    dx = 0.01;
    xf(1).value = 20;     % xf > x2 (always)
%     xf(1).value = 8;     % xf > x2 (always)
    xf(1).descriptor = 'Initial position of traffic state F';
    
    
%     x_cv2 = logspace(log10(x1(1).value),log10(xf(1).value), 200);
%     x_cv2 = linspace(x1(1).value, xf(1).value, 50);
    x_cv2 = linspace(x1(1).value, 2.1*x1(1).value, 250);
    
    
    %% Analytical calculation of horizons
    % Parameters
    pa = 1/ve;
    pb = (vf-ve)/(ve*(vah + vf));
    pc = (1/vh) * ((vah+vh) / (vah+vf));
    pd = 1/vf;
    pl = (1/vef) * ((vae+vef) / (vf+vae));
    pm = (1/vef) * ((vf-vef) / (vf+vae));
    pn = (1/vef) * ((vf - vef) / (vah + vf));
%     p_beta = qa*p_d - qi*p_c - qg*p_d;
    
    %Quadratic coefficients
    %D0 = 0.5*(qa-qi)*(vah/(vf*vfh))*((vah + vfh) / (vah + vf))*(xf(1).value)^2;
    D0(count_param) = 0.5*(qa-qi)*(1/(vf*vfh))*((vah + vfh) / (vah + vf))*(xf(1).value)^2;
    
    PA = 0.5*(qa-qg)*pl*pd*(xf(1).value)^2;
    PB = 0.5*(qi*pc*(pc-pm) - qg*pa*(pa-pm));
    PC = (0.5*qg*pb*pn - 0.5*qg*pb^2)*(x1(1).value)^2;
    
    PE = (qg*pa*pb - 0.5*(qg*pa - qi*pc)*pn - 0.5*qg*pb*pm)*x1(1).value;
    PF = (0.5*(qa-qg)*pd*pm + 0.5*(pl+pd)*(qg*pa-qi*pc))*xf(1).value;
    
    PG = (-qg*pb*pd - 0.5*(qa-qg)*pd*pn - 0.5*qg*pb*(pl-pd))*xf(1).value*x1(1).value;
    
    horz = roots([PB, PE + PF, PA+PC+PG-D0(count_param)]);
    horizon(:, count_param) = horz
%     horizon = roots([PB, PG*x1(1).value + PF, PC*(x1(1).value).^2 + PE*x1(1).value + PA]);
%     horizon(:, count_param) = horizon - x1(1).value;     
        
    if horizon(1, count_param) <= horizon(2, count_param)
        if horizon(1, count_param) > 0
            horizon_null(count_param) = horizon(1, count_param);
        else
            horizon_null(count_param) = horizon(2, count_param);
        end
%         horizon_event(count_param) = horizon(1, count_param);
    else
        if horizon(2, count_param) > 0
            horizon_null(count_param) = horizon(2, count_param);
        else
            horizon_null(count_param) = horizon(1, count_param);
        end
%         horizon_event(count_param) = horizon(2, count_param);
    end    
    
    horizon_event(count_param) = x1(1).value/(1+(vah/vh)*((vh-ve)/(vf-ve))) - x1(1).value;
    
    
    if flag_plot_dela == 1
        h_fig20 = figure(20);
        set(h_fig20, 'Position', [10.2*figProp figProp 4*figProp 3*figProp]);
        plot_x = linspace(x1(1).value, 1.2*x1(1).value, 50);
        %     plot_dela = PB.*plot_x.^2 + (PG.*x1(1).value + PF).*plot_x + PC.*x1(1).value.^2 + PE.*x1(1).value + PA;
        plot_dela = PB.*plot_x.^2 + (PE + PF).*plot_x + (PA + PC + PG - D0(count_param));
        plot(plot_x, plot_dela, 'Linewidth', 2);
        grid on;
    end
    
    
    % LOOP 2: Distance of CV2
    while count_delay < length(x_cv2)
        x2(1).value = x_cv2(count_delay);
        
        
        AR1(count_param, count_delay) = 0.5*(qa-qi)*pc^2*(x2(1).value)^2;
        AR2(count_param, count_delay) = (0.5*qa*(pa^2-pc^2) - qi*(pa*pc - pc^2))*(x2(1).value)^2 ...
            + 0.5*qa*pb^2*(x1(1).value)^2 + (qi*pb*pc - qa*pa*pb)*x1(1).value*x2(1).value;
        
        
        %     0.5*qa*pb^2*(x1(1).value)^2 + (0.5*qa*(pa^2-pc^2) ...
        %         - qi*(pa*pc - pc^2))*(x2(1).value)^2 + (qi*pb*pc-qa*pa*pb)*x1(1).value*x2(1).value;
        AR3(count_param, count_delay) = 0.5*(qa-qg)*pd^2*(xf(1).value)^2 ...
            - 0.5*(qa+qg)*pb^2*(x1(1).value)^2 + (-0.5*qa*pa^2 - 0.5*qg*pa^2 ...
            + qi*pa*pc)*(x2(1).value)^2 + (qa*pa*pb - qi*pb*pc + qg*pa*pb)*x1(1).value*x2(1).value ...
            + (-qi*pc*pd + qg*pa*pd)*xf(1).value*x2(1).value ...
            - qg*pb*pd*x1(1).value*xf(1).value;
        
        AR4(count_param, count_delay) = 0.5*((qa-qg)*pd*xf(1).value + (qg*pa-qi*pc)*x2(1).value ...
            -qg*pb*x1(1).value)*((pl-pd)*xf(1).value + pm*x2(1).value - pn*x1(1).value);
        
        if(x2(1).value < x1(1).value)
            error('Error - x2 must always be greater than x1');
        end
        if(xf(1).value < x2(1).value)
            error('Error - xf must always be greater than x2');
        end
              
        % Virtual arrival time calculations
        
        % Vertex 0-1
        T0(1).value = 0;
        x0(1).value = 0;
        vertex0(1).T = T0(1).value;
        vertex0(1).x = x0(1).value;
        
        % Vertex 0-2
        T0(2).value = xf(1).value/(vah + vf);
        T0(2).label = 'T_ahf';
        T0(2).descriptor = 'Time of intersection of A, H and F';
        x0(2).value = vah*T0(2).value;
        x0(2).descriptor = 'Position of first vehicle in state F (without CVs) when it joins the queue';
        vertex0(2).T = T0(2).value;
        vertex0(2).x = x0(2).value;
        
        % Vertex 0-3
        T0(3).value = T0(2).value + x0(2).value/vfh;
        T0(3).label = 'T_vanish';
        T0(3).descriptor = 'Time at which bottleneck vanishes';
        x0(3).value = 0;
        x0(3).descriptor = 'Bottlneck vanishes';
        vertex0(3).T = T0(3).value;
        vertex0(3).x = x0(3).value;
        
        %Vertex 0-4-arrival
        T0(4).value = xf(1).value/vf;
        T0(4).descriptor = 'Virtual arrival time for state F at bottleneck';
        x0(4).value = 0;
        x0(4).descriptor = 'Virtual arrival time for state F at bottleneck';
        vertex0(4).T = T0(4).value;
        vertex0(4).x = x0(4).value;
        

        
        % Vertex 2
        T(2).value = x1(1).value/(vf+vah);
        T(2).label = 'T_j^1';
        T(2).descriptor = 'Time when CV1 joins the queue';
        
        x1(2).value = vah*T(2).value;
        x1(2).descriptor = 'CV1 position when it joins the queue';
        vertex(2).x = x1(2).value;
        vertex(2).T = T(2).value;
        vertex(2).descriptor = 'CV1 joins queue';
        
        % Vertex 3
        T(3).value = x2(1).value/(vf+vah);
        T(3).label = 'T_j^2-';
        T(3).descriptor = 'Vehicle preceding CV2 joins the queue';
        
        x2(2).value = vah*T(3).value;
        x2(2).descriptor = 'Position of vehicle just downstream of CV2 when it joins the queue';
        vertex(3).x = x2(2).value;
        vertex(3).T = T(3).value;
        vertex(3).descriptor = 'Vehicle preceding CV2 joins the queue';
        
        % Vertex 4
        x2(3).value = x2(1).value - vf*T(2).value;
        x2(3).descriptor = 'Position of CV2 when it receives signal of CV1 joining the queue';
        vertex(4).x = x2(3).value;
        vertex(4).T = T(2).value;
        vertex(4).descriptor = 'CV2 receives CV1 signal';
        
        % Vertex 5
        T(4).value = ((vf-vh)*T(3).value - (vf-ve)*T(2).value)/(ve-vh);
        T(4).label = 'T_j^2';
        T(4).descriptor = 'Time when CV2 joins the queue';
        
        x2(4).value = x2(3).value - ve*(T(4).value - T(2).value);
        x2(4).descriptor = 'Position of CV2 when it joins the queue';
        vertex(5).x = x2(4).value;
        vertex(5).T = T(4).value;
        vertex(5).descriptor = 'CV2 joins the queue';
        
        % Case C - Vertex 11
        T(11).value = T(3).value + x2(2).value/vh;
        T(11).label = 'T_j^2-(alternate)';
        T(11).descriptor = 'Time when open zone begins';
        
        x2(7).value = 0;
        x2(7).descriptor = 'Bottleneck position';
        vertex(11).x = x2(7).value;
        vertex(11).T = T(11).value;
        vertex(11).descriptor = 'Vehicle before CV2 reaches bottleneck';
        
        % Case C - Vertex 14
        T(13).value = T(2).value + x2(3).value/ve;
        T(13).label = 'T_7''';
        T(13).descriptor = 'Time when CV2 reaches bottlneck (if it is farther upstream)';
        
        x2(8).value = 0;
        x2(8).descriptor = 'CV2 position at bottleneck';
        vertex(14).x = x2(8).value;
        vertex(14).T = T(13).value;
        vertex(14).descriptor = 'CV2 reaches the bottleneck (from farther upstream)';
        
        %Vertex 6
        T(5).value = T(2).value + x1(2).value/vh;
        T(5).label = 'T_l^1';
        T(5).descriptor = 'Time when CV1 leaves the queue/bottleneck';
        
        x1(3).value = 0;
        x1(3).descriptor = 'CV1 reaches bottleneck/leaves queue';
        vertex(6).x = x1(3).value;
        vertex(6).T = T(5).value;
        vertex(6).descriptor = 'CV1 reaches bottleneck';
        
        %Vertex 7
        T(6).value = T(4).value + x2(4).value/vh;
        T(6).label = 'T_l^2';
        T(6).descriptor = 'Time when CV2 leaves the queue/bottleneck';
        
        x2(5).value = 0;
        x2(5).descriptor = 'CV2 reaches bottleneck/leaves queue';
        vertex(7).x = x2(5).value;
        vertex(7).T = T(6).value;
        vertex(7).descriptor = 'CV2 reaches bottleneck';
        
        % Vertex 8
        T(7).value = ((xf(1).value - x2(1).value) + (vf + vae)*T(2).value)/(vf + vae);
        T(7).label = 'T7';
        T(7).descriptor = 'Time when State F hits state E';
        
        xf(2).value = xf(1).value - vf*T(7).value;
        xf(2).descriptor = 'Position of State F when it hits state E';
        vertex(8).x = xf(2).value;
        vertex(8).T = T(7).value;
        vertex(8).descriptor = 'State F hits state E';
        
        %Vertex 9
        T(8).value = (ve*(T(4).value - T(2).value) + w*T(4).value +...
            vae*(T(7).value - T(2).value) + vef*T(7).value)/(vef + w);
        T(8).label = 'T8';
        T(8).descriptor = 'Time when state F meets state H';
        
        xf(3).value = xf(2).value - vef*(T(8).value - T(7).value);
        xf(3).descriptor = 'Position of state F when it hits State H';
        vertex(9).x = xf(3).value;
        vertex(9).T = T(8).value;
        vertex(9).descriptor = 'State F meets state H';
        
        
        % Vertex 12
        T(10).value = ((vae + ve)*T(2).value - (ve + w)*T(4).value)/(vae - w);
        T(10).label = 'Tf-h';
        T(10).descriptor = 'State A hits state H, before F hits E';
        x2(6).value = x2(3).value + vae*(T(10).value - T(2).value);
        x2(6).descriptor = 'State F is very close to CV2';
        vertex(12).x = x2(6).value;
        vertex(12).T = T(10).value;
        vertex(12).descriptor = T(10).descriptor;
        
        
        
        
        %Vertex 10
        T(9).value = T(8).value + xf(3).value/vfh;
        T(9).label = 'Tf';
        T(9).descriptor = 'Time when state F reaches bottleneck';
        
        xf(4).value = 0;
        xf(4).descriptor = 'Position of state F when it hits the bottleneck point';
        vertex(10).x = xf(4).value;
        vertex(10).T = T(9).value;
        vertex(10).descriptor = 'State F hits bottleneck';
        
        %Vertex 13
        T(12).value = T(7).value + xf(2).value/vef;
        T(12).label = 'Tf''';
        T(12).descriptor = 'Time when state F reaches bottleneck (from far away)';
        
        xf(5).value = 0;
        xf(5).descriptor = 'Position of state F when it hits the bottleneck point';
        vertex(13).x = xf(5).value;
        vertex(13).T = T(12).value;
        vertex(13).descriptor = 'State F hits bottleneck (from far away)';
        
        
        %% Plotting
        c = [0.2 0.2 0.2];
        
        if flag_plot == 1
            
            h_fig1 = figure(1);
%             set(h_fig1, 'Position', [2*figProp 4.8*figProp 4*figProp 3*figProp]);
            set(h_fig1, 'Position', [2*figProp 2.8*figProp 5*figProp 4*figProp]);
            clf(h_fig1,'reset');
            ylim([-10 0])
            
            grid on;
            hold on;
            h = xlabel('Time (hr)');
            h = ylabel('Distance (km)');
            
            
            % Color the areas and put into background
            % Area in front of CV1
            h = fill([horzcat(T([1,3]).value),0], [-1.*[x2(1:2).value], 0],...
                [0.4 1 0.4],'EdgeColor','none');
            alpha(h,0.6);
            % Area between CV1 and CV2
            h = fill([0, vertex(4).T, vertex(8).T, 0], [-x2(1).value, -vertex(4).x,...
                -vertex(8).x, -xf(1).value], [0.4 1 0.4], 'EdgeColor','none');
            alpha(h,0.6);
            % Queued area
            h = fill([0, vertex(3).T, vertex(7).T], [0, -vertex(3).x, -vertex(7).x],...
                [1 0.4 0.4],'EdgeColor','none');
            alpha(h,0.6);
            
            if vertex(5).x > 0
                
                if vertex(9).x <= vertex(8).x
                    % Polygon 5-7-10-9-5 - State H
                    h = fill([vertex(5).T, vertex(7).T, vertex(10).T, vertex(9).T],...
                        [-vertex(5).x, -vertex(7).x, -vertex(10).x, -vertex(9).x], ...
                        [1 0.4 0.4], 'EdgeColor','none');
                    alpha(h,0.6);
                    
                    % Polygon 4-5-9-8-4 - State E
                    h = fill([vertex(4).T, vertex(5).T, vertex(9).T, vertex(8).T],...
                        [-vertex(4).x, -vertex(5).x, -vertex(9).x, -vertex(8).x], ...
                        [1 1 0.4], 'EdgeColor','none');
                    alpha(h,0.6);
                else
                    % Polygon 5-7-10-8-12-5 - State H
                    h = fill([vertex(5).T, vertex(7).T, vertex(10).T, vertex0(2).T, vertex(12).T],...
                        [-vertex(5).x, -vertex(7).x, -vertex(10).x, -vertex0(2).x, -vertex(12).x], ...
                        [1 0.4 0.4], 'EdgeColor','none');
                    alpha(h,0.6);
                    
                    % Polygon 4-5-12-4 - State E
                    h = fill([vertex(4).T, vertex(5).T, vertex(12).T],...
                        [-vertex(4).x, -vertex(5).x, -vertex(12).x], ...
                        [1 1 0.4], 'EdgeColor','none');
                    alpha(h,0.6);
                end
            else
                h = fill([vertex(4).T, vertex(14).T, vertex(13).T, vertex(8).T],...
                    [-vertex(4).x, -vertex(14).x, -vertex(13).x, -vertex(8).x], ...
                    [1 1 0.4], 'EdgeColor','none');
                alpha(h,0.6);
            end
            
            
            
            % Plotting x-t diagram lines in absence of CVs
            if flag_plot_noCV == 1
                gray = [0.5 0.5 0.5];
                line([vertex0(1).T, vertex0(2).T], [-vertex0(1).x, -vertex0(2).x], 'Linewidth', 2, 'Linestyle','--','Color', gray);
                line([vertex0(2).T, vertex0(3).T], [-vertex0(2).x, -vertex0(3).x], 'Linewidth', 2, 'Linestyle','--','Color', gray);
                line([vertex0(2).T, vertex0(4).T], [-vertex0(2).x, -vertex0(4).x], 'Linewidth', 2, 'Linestyle',':','Color', gray);
            end
            
            line([T(1).value, T(2).value], [0, -x1(2).value], 'Linewidth', 2, 'Color', c);
            % Line joining vertices 1 and 2
            % Line joining vertices x1 and 2
            line([T(1).value, T(2).value], [-x1(1).value, -x1(2).value],'Linewidth', 2, 'Color', 'red');
            
            % Line joining vertices 1 and 3
            line([T(1).value, T(3).value], [0, -x2(2).value], 'Linewidth', 2, 'Color', c);
            
            % Line joining initial x2 to vertex 4, and onwards to 3
            line([T(1).value, vertex(4).T], [-x2(1).value, -vertex(4).x],'Linewidth', 2, 'Color', 'red');
            line([vertex(4).T, vertex(3).T], [-vertex(4).x, -vertex(3).x], 'Linewidth', 2, 'Color', c);
            
            % Line joining vertices 2 and 6
            line([vertex(2).T, vertex(6).T], [-vertex(2).x, -vertex(6).x], 'Linewidth', 2, 'Color', 'red');
            
            % Line joining state F at time t = 0 to vertex 8,
            line([0, vertex(8).T], [-xf(1).value, -vertex(8).x], 'Linewidth', 2, 'Color', c);
            
            
            if vertex(5).x > 0
                % Line joining vertices 4 and 5
                line([vertex(4).T, vertex(5).T], [-vertex(4).x, -vertex(5).x], 'Linewidth', 2, 'Color', 'red');
                
                % Line joining vertices 3 to 5, and 5 to 7
                line([vertex(3).T, vertex(5).T], [-vertex(3).x, -vertex(5).x], 'Linewidth', 2, 'Color', c);
                line([vertex(5).T, vertex(7).T], [-vertex(5).x, -vertex(7).x], 'Linewidth', 2, 'Color', 'red');
                
                if vertex(9).x < vertex(8).x
                    % Line joining vertex 8 to 9
                    line([vertex(8).T, vertex(9).T], [-vertex(8).x, -vertex(9).x], 'Linewidth', 2, 'Color', c);
                    % Line joining vertex 4 to vertex 8
                    line([vertex(4).T, vertex(8).T], [-vertex(4).x, -vertex(8).x], 'Linewidth', 2, 'Color', c);
                    % Line joining vertex 4 to vertex 8
                    line([vertex(5).T, vertex(9).T], [-vertex(5).x, -vertex(9).x], 'Linewidth', 2, 'Color', c);
                    % Line joining vertex 9 to 10
                    line([vertex(9).T, vertex(10).T], [-vertex(9).x, -vertex(10).x], 'Linewidth', 2, 'Color', c);
                else
                    % Line joining vertex 4 to vertex 12
                    line([vertex(4).T, vertex(12).T], [-vertex(4).x, -vertex(12).x], 'Linewidth', 2, 'Color', c);
                    % Line joining vertex 5 to vertex 12
                    line([vertex(5).T, vertex(12).T], [-vertex(5).x, -vertex(12).x], 'Linewidth', 2, 'Color', c);
                    % Line joining vertex 12 to vertex0(2)
                    line([vertex(12).T, vertex0(2).T], [-vertex(12).x, -vertex0(2).x], 'Linewidth', 2, 'Color', c);
                    % Line joining original F-H interface
                    line([vertex0(2).T, vertex0(3).T], [-vertex0(2).x, -vertex0(3).x], 'Linewidth', 2, 'Color', c);
                end
            else
                % Line joining vertex 3 to 11
                line([vertex(3).T, vertex(11).T], [-vertex(3).x, -vertex(11).x], 'Linewidth', 2, 'Color', c);
                % Line joining vertex 4 to 14
                line([vertex(4).T, vertex(14).T], [-vertex(4).x, -vertex(14).x], 'Linewidth', 2, 'Color', 'red');
                % Line joining vertex 4 to 8
                line([vertex(4).T, vertex(8).T], [-vertex(4).x, -vertex(8).x], 'Linewidth', 2, 'Color', c);
                % Line joining vertex 8 to 13 (10')
                line([vertex(8).T, vertex(13).T], [-vertex(8).x, -vertex(13).x], 'Linewidth', 2, 'Color', c);
            end
        end
        
        
        %% Delay curve at bottleneck
        if vertex(5).x > 0
            delay = polyarea([0, vertex0(4).T, vertex0(3).T], [0, ...
                qa*vertex0(4).T, qh*vertex(10).T]);
        else
            delay = polyarea([0, vertex0(4).T, vertex(13).T, vertex(14).T, ...
                vertex(11).T], [0, qa*vertex0(4).T, qh*vertex(11).T + ...
                qe*(vertex(13).T - vertex(14).T), qh*vertex(11).T,...
                qh*vertex(11).T]);
        end
        
        if flag_plot_N_t == 1
            
            h_fig5 = figure(5);
            %         set(h_fig5, 'Position', [2*figProp figProp 4*figProp 3*figProp]);
            set(h_fig5, 'Position', [2*figProp figProp 5*figProp 4*figProp]);
            clf(h_fig5,'reset');
            
            h = xlabel('Arrival/Departure time at bottleneck (hr)');
            h = ylabel('Number of vehicles');
            h = title('N-t curve at bottleneck');
            grid on;
            hold on;
            % Virtual arrival curve
            line([0, vertex0(4).T], [0, qa*vertex0(4).T], 'Linewidth', 2, 'Color', c)
            line([vertex0(4).T, vertex0(3).T], [qa*vertex0(4).T, qa*vertex0(4).T + qf*(vertex0(3).T-vertex0(4).T)], 'Linewidth', 2, 'Color', c)
            
            % Departure curve with no CVs
            line([0, vertex(10).T], [0, qh*vertex(10).T], 'Linewidth', 2, ...
                'Linestyle', ':', 'Color', c)
            
            % Departure curve
            if vertex(5).x > 0
                line([0, vertex(10).T], [0, qh*vertex(10).T], 'Linewidth', 2, ...
                    'Color', c)
                line([0, vertex(10).T], [0, qh*vertex(10).T], 'Linewidth', 2, ...
                    'Linestyle', ':', 'Color', c)
                
                h = fill([0, vertex0(4).T, vertex0(3).T], [0, qa*vertex0(4).T, ...
                    qh*vertex(10).T], [0.4 0.4 1], 'EdgeColor','none');
                alpha(h,0.3);
                
                delay = polyarea([0, vertex0(4).T, vertex0(3).T], [0, ...
                    qa*vertex0(4).T, qh*vertex(10).T]);
                
            else
                line([0, vertex(11).T], [0, qh*vertex(11).T], 'Linewidth', 2, ...
                    'Color', c)
                line([vertex(11).T, vertex(14).T], [qh*vertex(11).T, qh*vertex(11).T], 'Linewidth', 2, ...
                    'Color', c)
                line([vertex(14).T, vertex(13).T], [qh*vertex(11).T, qh*vertex(11).T...
                    + qe*(vertex(13).T - vertex(14).T)], 'Linewidth', 2,'Color', c)
                
                
                h = fill([0, vertex0(4).T, vertex(13).T, vertex(14).T, vertex(11).T], [0, qa*vertex0(4).T, ...
                    qh*vertex(11).T + qe*(vertex(13).T - vertex(14).T), qh*vertex(11).T, qh*vertex(11).T],...
                    [0.4 0.4 1], 'EdgeColor','none');
                alpha(h,0.3);
                
                delay = polyarea([0, vertex0(4).T, vertex(13).T, vertex(14).T, ...
                    vertex(11).T], [0, qa*vertex0(4).T, qh*vertex(11).T + ...
                    qe*(vertex(13).T - vertex(14).T), qh*vertex(11).T,...
                    qh*vertex(11).T]);
            end
        end
               
        % Save delay data
        delay_data(count_param, count_delay) = delay;
        loc(count_delay) = x2(1).value - x1(1).value;
        count_delay = count_delay + 1;
                
%        pause
        
%         if vertex(13).T > vertex0(3).T
%             break
%         end       
    end
    delay_data(count_param,:) = delay_data(count_param,:) - delay_data(count_param,1);    
    
    if flag_plot_delay_profile == 1
        % Figure for total delay data
        h_fig10 = figure(10);
        set(h_fig10, 'Position', [6.1*figProp 4.8*figProp 4*figProp 3*figProp]);
        plot(loc, delay_data, 'Linewidth',2);
        h = xlabel('Distance between CV1 and CV2 (km)');
        h = ylabel('Total delay (veh-hrs)');
        h = title('Total delay at bottleneck');
        grid on;
        hold on;
    end
    
    count_param = count_param + 1
 
end


%% Plotting influential subspaces

h_fig15 = figure(15);
set(h_fig15, 'Position', [6.1*figProp figProp 4*figProp 3*figProp]);
% For changing ka
% delay_data((delay_data == 0)) = NaN;

delay_data(abs(delay_data) < 5e-5) = 0;
levels = linspace(min(min(delay_data)), max(max(delay_data)),21);

% levels = logspace(-2,2,100);
% levels = linspace(-2,15,20);
% levels = [-0.3, 0, 1, 3, 5, 7, 9, 11, 13];
levels = [-250, -50, -20, -10, 0, 10, 20, 30, 40, 50, 100];

% contourf(loc, kf_data./kc, delay_data, levels, 'Showtext', 'on');
% levels = [-150, -25, -10, -5, 0, 5, 10, 25, 50,75];
contourf(loc, ka_data./kc, delay_data, levels, 'Showtext', 'on');
% levels = [-150, -25, -10, -5, -2, 0, 5, 10, 25, 50];
% contourf(loc, kg_data./ka, delay_data, levels, 'Showtext', 'on');
h = xlabel('Distance between CV1 and CV2 (km)');
h = ylabel('Upstream traffic density (k_F/k_C)');

hold on;
% plot(abs(horizon_event), ka_data./kc, 'r-', 'Linewidth',2);
% plot(abs(horizon_null-x1(1).value), ka_data./kc, 'b-', 'Linewidth',2);
plot(abs(horizon_event), kf_data./kc, 'r-', 'Linewidth',2);
plot(abs(horizon_null-x1(1).value), kf_data./kc, 'b-', 'Linewidth',2);
% plot(abs(horizon_event), kg_data./ka, 'r-', 'Linewidth',2);
% plot(abs(horizon_null-x1(1).value), kg_data./ka, 'b-', 'Linewidth',2);
grid on;
title(['x_1 = ', num2str(x1(1).value), ' km,  and x_f = ',num2str(xf(1).value), ' km']);