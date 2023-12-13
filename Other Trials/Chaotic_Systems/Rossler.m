function y = Rossler(t,x)

% parameters

a = 0.25; b = 0.2; c = 1;

% equations

y = [
     -x(2)-x(3)+sin(t);            % dx/dt
     x(1)+a*x(2);           % dy/dt
     b+x(3)*(x(1)-c);       % dz/dt
] ; 
