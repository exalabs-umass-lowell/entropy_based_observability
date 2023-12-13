clc;
clear;
%% Initial conditions and parameters
x0 = -10;
x1 = -8;


% Establish the van der pol model
t  = 0:0.1:500;   % time scale
[t,xi] = ode45( 'vd', t, [x0 x1]); 

x = interp1(xi,t,'linear');

N = size(x,1);
S = zeros(N, N);

for i = 1:N,
    S(:,i) = abs( repmat( x(i), N, 1 ) - x(:) );
end

imagesc(t, t, S)
axis square