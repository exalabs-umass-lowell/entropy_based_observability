
function y=Neuron_Three_HR(t,M)

%Initial conditiosn and parameter values

a = 1;
b = 3.0;
c = 1.0;
d = 5.0;
s = 4.0;
r = 0.006;
x0 = -1.6;
Epsilon = 0.5;
B3 = 500;
B2 = 50;
B1 = 10;
x1 = M(1);
y1 = M(2);
z1 = M(3);

x2 = M(4);
y2 = M(5);
z2 = M(6);

x3 = M(7);
y3 = M(8);
z3 = M(9);

I_ext = 3;
Excit = 10;
y = [B1*(y1-a*x1^3+b*x1^2-z1+I_ext+Epsilon*(x2-x1));
    B1*(c-d*x1^2-y1);
    (s*(x1-x0)-z1)+Excit;
    B2*(y2-a*x2^3+b*x2^2-z2+I_ext+Epsilon*((x1-x2)+(x3-x2)));
    B2*(c-d*x2^2-y2);
    (s*(x2-x0)-z2);
    B3*(y3-a*x3^3+b*x3^2-z3+I_ext+Epsilon*(x2-x3));
    B3*(c-d*x3^2-y3);
    (s*(x3-x0)-z3)];
    