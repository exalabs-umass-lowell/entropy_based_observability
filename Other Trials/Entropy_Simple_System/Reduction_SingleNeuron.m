% Here we use order reduction method to reduce the system to two
% independent first order system




function y = Reduction_SingleNeuron(t,x)

% Parameters
a = -1;
b = -2;
k = 0.006;
Epsilon = 0.5; 
c = -Epsilon;

% Simplification by replacing with other coefficients 
d = 4-x(1);
p = -b/3/a;
q = p^3+(b*c-3*a*d)/6/a^2;
r = c/3/a;

d2 = 4-x(2);
p2 = -b/3/a;
q2 = p^3+(b*c-3*a*d2)/6/a^2;
r = c/3/a;

s = (q+(q^2+(r-p^2)^3)^(1/2))^(1/3) + (q-(q^2+(r-p^2)^3)^(1/2))^(1/3) + p;

% Simplification by replacing with other coefficients 

s2 = (q+(q2^2+(r-p2^2)^3)^(1/2))^(1/3) + (q2-(q2^2+(r-p2^2)^3)^(1/2))^(1/3) + p2;
y = [k*(4.0*(s+1.60)-x(1))+5;
    k*(4.0*(s2+1.60)-x(2))];

