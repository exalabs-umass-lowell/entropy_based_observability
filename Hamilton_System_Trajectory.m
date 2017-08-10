% Simulation under different epsilon to see the trajectories of Hamiltonian
% system


x = -pi:0.02:pi;
epsilon = [0.017, 0.027, 0.037];
for i = 1:3
    figure;
    for t = 0:10
        p = sqrt(2*epsilon(i)*(cos(x)+cos(x-t)));
        
        plot(x,p,x,-p);
        hold on;
    end
end