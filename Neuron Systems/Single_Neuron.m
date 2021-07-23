function f = Single_Neuron(t,X)

a = 1;
b = 3.0;
c = 1.0;
d = 5.0;
s = 4.0;
x0 = -1.6;
I_ext = 3;


x1 = X(1);
y1 = X(2);
z1 = X(3);
f = [(y1-a*x1^3+b*x1^2-z1+I_ext);
     (c-d*x1^2-y1);
     0.01*(s*(x1-x0)-z1)+5];