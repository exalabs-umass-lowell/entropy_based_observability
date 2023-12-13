function y = HyperRossler_SingularPerturb(t,x)

a= 0.25;
b = 3;
k =  0.5;

y = [-x(2) - x(3);
    x(1) + a*x(2)+10*x(3);
    b + x(1)*x(3)];