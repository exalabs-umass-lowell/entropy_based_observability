function f = Reduced_Single_Neuron(t,X)

a = 1;
b = 3.0;
c = 1.0;
d = 5.0;
s = 4.0;
x0 = -1.6;
I_ext = 3;

x1 = X(1);
z1 = X(2);
f = [ ((s*(x1-x0)-z1))/(0.1+2*(b-d)*x1-3*a*x1^2);
      (s*(x1-x0)-z1)+5];