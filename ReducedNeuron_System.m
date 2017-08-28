function f=ReducedNeuron_System(t,M)

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
z1 = M(2);

x2 = M(3);
z2 = M(4);

I_ext = 3.0;

f = [r*(s*(x1-x0)-z1)/((b-d)*x1*2-3*a*x1^2-Epsilon+Epsilon*r*(s*(x2-x0)-z2)/((b-d)*x2*2-3*a*x2^2-Epsilon));
    r*(s*(x1-x0)-z1)+5;
    r*(s*(x2-x0)-z2)/((b-d)*x2*2-3*a*x2^2-Epsilon+Epsilon*r*(s*(x1-x0)-z1)/((b-d)*x1*2-3*a*x1^2-Epsilon));
    r*(s*(x2-x0)-z2)];
    