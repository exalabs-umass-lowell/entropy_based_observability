function [ robot_heading_new ] = swarm_fn(robot_pos_old, robot_heading_old, param)
% Function to output new heading of robots based on old positions, old
% robot headings and parameter values

% Define the parameters
N = param.N;
dT = param.dT;
r_repulsion = param.r_repulsion;                % Zone of repulsion radius
r_orientation = param.r_orientation;            % Zone of orientation radius
r_attraction = param.r_attraction;              % Zone of attraction radius
perception_field = param.perception_field;      % Field of view of the robot [0, 360] (degrees)
blind_zone = (360 - perception_field)/2;
theta_dot = param.theta_dot;                    % Turning rate (deg/s)
speed = param.speed;                            % Robot speed (m/s) - All robots possess same velocity


% Identify the distances to neigboring robots
robot_distance = pdist2(robot_pos_old,robot_pos_old);
robot_direction = zeros(N,2);

% NOTE: There is definitely a better way of doing this
for(i = 1:1:N)
    for(j = 1:1:N)
        % Identify the heading to the neighboring robots
        angle(i,j) = atan2((robot_pos_old(i,2) - robot_pos_old(j,2)),(robot_pos_old(i,1) - robot_pos_old(j,1)));
    end
    robot_direction(i,1) = cosd(robot_heading_old(i));
    robot_direction(i,2) = sind(robot_heading_old(i));
end
angle = 180 + angle.*(180/pi);    % Convert radians to degrees
% Identify the robots that are in the field of perception
[vis_row, vis_col] = find(angle - blind_zone > 0 & angle < 360 - blind_zone);
temp_row = vis_row(vis_row - vis_col ~= 0);
temp_col = vis_col(vis_row - vis_col ~= 0);
vis_row = temp_row;
vis_col = temp_col;
dr = zeros(N,2);
di = zeros(N,2);
% If neighbors present in zone of repulsion (ZoR), use avoid algorithm
[rep_row, rep_col] = find(robot_distance > 0 & robot_distance < r_repulsion);
for i=1:size(rep_row)
    r_ij = (robot_pos_old(rep_row(i),:) - robot_pos_old(rep_col(i),:))/norm(robot_pos_old(rep_row(i),:) - robot_pos_old(rep_col(i),:));
    dr(rep_col(i),:);
    dr(rep_col(i),:) = dr(rep_col(i),:) + r_ij;%/norm(r_ij);
    norm(r_ij);
    %end
end

% If no neighbors in zone of repulsion, then
    % If neighbors only present in zone of orientation (ZoO), use rule d_o
[or_row, or_col] = find(robot_distance>r_repulsion & robot_distance<r_orientation);
do = zeros(N,2);
    for i=1:size(or_col) 
        if dr(or_col(i),1) == 0 %&& norm(robot_orientation(or_col(i),:) - robot_orientation(or_row(i),:)) ~= 0
            v_ij = (robot_direction(or_row(i),:))/norm(robot_direction(or_row(i),:));
            do(or_col(i),:) =  do(or_col(i),:) + v_ij;%/norm(v_ij);
        end
    end

[at_row, at_col] = find(robot_distance<r_attraction & robot_distance>r_orientation); 
da = zeros(N,2);
    for i=1:size(at_col) 
        if dr(at_col(i),1) == 0
            i;
            r_ij = (robot_pos_old(at_row(i),:) - robot_pos_old(at_col(i),:))/norm(robot_pos_old(at_row(i),:) - robot_pos_old(at_col(i),:));
            da(at_col(i),:) =  da(at_col(i),:) + r_ij;%/norm(r_ij);
        end
    end
do;
dr = -dr;
da;
    % If neighbors only present in zone of attraction (ZoA), use rule d_a
    for i=1:N
        if da(i,1) ~=0 && do(i,1) ~=0
            di(i,:) = 0.5*(da(i,:)+do(i,:));
            %di(i,:) = di(i,:)/norm(di(i,:));
        elseif da(i,1) ==0 && do(i,1) ~=0
            di(i,:) = do(i,:);
            %di(i,:) = di(i,:)/norm(di(i,:));
        elseif da(i,1) ~=0 && do(i,1) ==0
            di(i,:) = da(i,:);
            %di(i,:) = di(i,:)/norm(di(i,:));
        %else
            %di(i,1) = cosd(robot_heading_old(i));
            %di(i,2) = sind(robot_heading_old(i));
        %end
        end
    di = di + dr;
    new_heading = zeros(N,1);
    end
    for i=1:N
        if di(i,:) ~=0
            new_heading(i,1) = atan2d(di(i,2),di(i,1)); %-10 + (10)*rand(1);
        end
    end
    di;
    new_heading;
    robot_pos_old;
    da;
    % If neighbors in both ZoO and ZoA, use rule 0.5*(d_o + d_a)

    % If social forces create zero vector, or no individuals are found,
    % use rule v (just move with previous velocity)

    % Add randomness via a Gaussian distribution to determine final robot heading


    % Use turning rate \theta to turn toward determined robot heading

%}

robot_heading_new = new_heading ; % for now!!!