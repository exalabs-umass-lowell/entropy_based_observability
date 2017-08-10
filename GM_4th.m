
% -------Brandon Yang

function f = GM_4th(x,n,C)

% This would output the function values of the IDM for N cars
% x(1), x(2),...x(N) are positions, and x(n+1), x(n+2),...x(2n) are
% velocities, x(2n+1), x(2n+2) etc are 1/x(1), 1/x(2)....

alpha = 0.7;
f = zeros(2*n,1);
v = zeros(n,1);
h = zeros(n,1);


%% System dynamics for the 1st vehicles

h(1) = C-(x(1)-x(n));           % Relative displacement
v(1) = x(2*n)-x(n+1);           % Relative velocity

f(1) = x(n+1);
f(n+1) = alpha*x(n+1)*v(1)/h(1);



%% System dynamics for the following vehicles

for i = 2:n
    hOld = x(i-1)-x(i);         % Relative displacement
    dv = x(n+i-1)-x(i+n);          % Relative velocity
    f(i) = x(n+i);
    f(n+i) = alpha*(dv)*x(n+i)/hOld;

end

end



