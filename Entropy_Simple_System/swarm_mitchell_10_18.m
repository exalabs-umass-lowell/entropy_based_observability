
%% Script to generate swarming behavior in robots in 2-dimensions
% Authors: Mitchell Scott and Kshitij Jerath
% Contact: mitchell.scott@wsu.edu (put your permanent email id here)


%%%%Script WITHOUT blindness and turning angle%%%%
%% Clear the workspace
clear all
clc
%clf reset
%% Define the model parameters
L  = 7.0;         % Cell size (m)
N = 100;         % Number of robots in swarm
param.N = N;    % Store as a parameter struct for easy passing
dT = 0.1;      % Simulation time increment (s)
param.dT = dT;  % Store as a parameter struct for easy passing
% density rho = N/L^2

%%%Some milling with L=7, N=100, r=1, o=16, a=20. Works with L from
%%%3-7 .5 DT == L=7, N=100, r=1, o=12, a =20

%%%Swarm: L=3, N=100, r=1, o=2, a=15

%%Parallel: L=5, N=100, r=1, o=14, a=15

% Define the radii for the interaction zones
param.r_repulsion = 1.0;
param.r_orientation = 9.75;
param.r_attraction = 20;
param.perception_field = 250;     % Field of view of the robot [0, 360] (degrees)
param.theta_dot = 10;             % Turning rate (deg/s)

param.speed = 1.8;                % robot speed (m/s) - All robots possess same velocity
speed = param.speed;
%% Initialize the system
num_iter = 1000;             % Number of iterations
sim_time = num_iter*dT;     % Total simulation time (s)

% Array size is number of (robots x space dimensions x iterations)
robot_pos = zeros(N,2,num_iter); % robot positions (all initialized to origin)

% Initialize the robot positions in the 2-D world according to a uniform
% distribution
robot_pos(:,1,1) = -L + (2*L)*rand(1,N);
robot_pos(:,2,1) = -L + (2*L)*rand(1,N);

% Initialize robot orientation in the 2-D world
robot_heading = zeros(N,1,num_iter);
robot_heading(:,1,1) = 360.*rand(N,1);    % Uniformly distributed orientations

% Initialize robot velocities
% Array size is number of (robots x space dimensions x iterations)
robot_vel = zeros(N,2,num_iter);
robot_vel(:,1,1) = speed.*cos(robot_heading(:,1,1));
robot_vel(:,2,1) = speed.*sin(robot_heading(:,1,1));

%print(1)
% Initialize direction vectors for zones
for iter = 2:1:num_iter
    clf
    robot_distance = pdist2(robot_pos(:,1,num_iter),robot_pos(:,1,num_iter));
    %[rep_row, rep_col] = find(robot_distance>0 & robot_distance<r_repulsion);    
    %[or_row, or_col] = find(robot_distance>r_repulsion & robot_distance<r_orientation);
    %[at_row, at_col] = find(robot_distance<r_attraction & robot_distance>r_orientation);
    % Identify new robot heading based on behavioral rules
    %di = swarm_fn_mit(robot_pos(:,:,iter-1), robot_heading(:,1,iter-1), param);
    
    robot_heading(:,1,iter) = swarm_fn_mit(robot_pos(:,:,iter-1), robot_heading(:,1,iter-1), param);
   
       
    
    % Include the dynamics of all robots - positions and velocities
    robot_pos(:,1,iter) = robot_pos(:,1,iter-1) + dT*robot_vel(:,1,iter-1);
    robot_pos(:,2,iter) = robot_pos(:,2,iter-1) + dT*robot_vel(:,2,iter-1);
    robot_vel(:,1,iter) = speed.*cosd(robot_heading(:,1,iter));
    robot_vel(:,2,iter) = speed.*sind(robot_heading(:,1,iter));
        
    % Need to add global metrics p_group, m_group and c_group
    
    
    % Plot the results
    quiver(robot_pos(:,1,iter), robot_pos(:,2,iter), robot_vel(:,1,iter), robot_vel(:,2,iter), 0.3, 'color',[0 0 0]);
    if iter < 10
        plot_until = iter - 1;
    else
        plot_unit = 10;
    end
    hold on
    if iter >1
        quiver(robot_pos(:,1,iter-plot_until:iter-1), robot_pos(:,2,iter-plot_until:iter-1), robot_vel(:,1,iter-plot_until:iter-1), robot_vel(:,2,iter-plot_until:iter-1),0.2, 'color',[0.5 0.5 1.0]);
    end
    xlabel('X position')
    ylabel('Y position')
    xlim([-30 30])
    ylim([-30 30])
    %scatter(robot_pos(:,1,iter), robot_pos(:,2,iter))
    %dT*iter
    pause(.001)

    
end
clf
pause(3)
%{
for iter = 2:4:num_iter
    clf
    
        quiver(robot_pos(:,1,iter), robot_pos(:,2,iter), robot_vel(:,1,iter), robot_vel(:,2,iter), 0.3, 'color',[0 0 0]);
    if iter < 20
        plot_until = iter - 1;
    else
        plot_unit = 20;
    end
    hold on
    if iter >1
        quiver(robot_pos(:,1,iter-plot_until:iter-1), robot_pos(:,2,iter-plot_until:iter-1), robot_vel(:,1,iter-plot_until:iter-1), robot_vel(:,2,iter-plot_until:iter-1),0.2, 'color',[0.5 0.5 1.0]);
    end
    
    xlabel('X position')
    ylabel('Y position')
    xlim([-30 30])
    ylim([-30 30])
    pause(0.1)
    iter
end
%}