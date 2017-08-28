% Models for two neuron system


function y=Neuron_System(t,M)

%Initial conditiosn and parameter values

a = 1;
b = 3.0;
c = 1.0;
d = 5.0;
s = 4.0;
r = 0.006;
x0 = -1.6;
Epsilon = 3;

x1 = M(1);
y1 = M(2);
z1 = M(3);

x2 = M(4);
y2 = M(5);
z2 = M(6);

I_ext = 3;

y = [y1-a*x1^3+b*x1^2-z1+I_ext+Epsilon*(x2-x1);
    c-d*x1^2-y1;
    r*(s*(x1-x0)-z1)+5;
    y2-a*x2^3+b*x2^2-z2+I_ext+Epsilon*(x1-x2);
    c-d*x2^2-y2;
    r*(s*(x2-x0)-z2)];
    