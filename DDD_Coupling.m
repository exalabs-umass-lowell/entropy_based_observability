% Referenced the paper: Effect of mixed coupling on relay-coupled R ?ossler and Lorenz oscillators

% The case for DDD

function y=DDD_Coupling(t,M)

x1 = M(1);
y1 = M(2);
z1 = M(3);

x2 = M(4);
y2 = M(5);
z2 = M(6);

x3 = M(7);
y3 = M(8);
z3 = M(9);

a = 0.1;
b = 0.1;
c = 14;
k = 0.1;

% g1 = k*(y2-x1);
% g2 = k*(y1-x2)+k*(y3-x2);
% g3 = k*(y2-x3);



g1 = k*(y2-x1);
g2 = k*(y1-x2)+k*(y3-x2);
g3 = k*(y2-x3);



y = [-y1-z1+g1+5;
     x1+y1;
     b+z1*(x1-c);

    -y2-z2+g2;
    x2+y2;
    b+z2*(x2-c);


    -y3-z3+g3;
    x3+y3;
    b+z3*(x3-c)];