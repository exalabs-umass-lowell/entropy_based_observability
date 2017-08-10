function y = HyperRossler(t,x)

a= 0.25;
b = 3;
c = 50;
d = 5;

y = [-x(2) - x(3);
    x(1) + a*x(2) + x(4);
    b + x(1)*x(3);
    -c*x(3)+d*x(4)];